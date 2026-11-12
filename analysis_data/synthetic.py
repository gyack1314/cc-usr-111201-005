"""合成对照（synthetic control）与空间 placebo 诊断。

本模块只做一件事：在一个平衡的周度厂区面板上，用上线前周拟合一个
非负、和为 1、无截距的 donor 权重，并给出上线后的差距曲线、留出周
选 lambda、以及把每个 donor 轮流当伪处理厂区时的空间 placebo 诊断。

处理后的数据只在权重冻结之后用于计算“观测值 - 合成值”。它不参与
标准化、不参与 lambda 选择、也不参与权重拟合。合成对照不能自动消除
未观察混杂；这里给出的是在该权重下的对照差距，不是因果证明。
"""
import math
from pathlib import Path

from .io import load

ACTIVE_TOL = 1e-12
DEFAULT_TIE_TOLERANCE = 1e-12
DEFAULT_ZERO_TOLERANCE = 1e-12
DEFAULT_SOLVER_TOLERANCE = 1e-8
DEFAULT_MAX_ITERATIONS = 20000

_POSITIVE_INFINITY = "+infinity"


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _number(value, name, *, nonnegative=False, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("%s 必须是有限数值" % name)
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("%s 必须是有限数值" % name)
    if nonnegative and value < 0:
        raise ValueError("%s 不能为负" % name)
    if positive and value <= 0:
        raise ValueError("%s 必须为正" % name)
    return value


def _mse(observed, weights, donor_series):
    """按 donor 序列和权重计算 observed 的均方误差。"""
    total = 0.0
    for t, target in enumerate(observed):
        synthetic = sum(weight * series[t] for weight, series in zip(weights, donor_series))
        residual = target - synthetic
        total += residual * residual
    return total / len(observed) if observed else float("nan")


def _fix_sum(weights):
    """把权重和校正到 1；校正量只会落在当前最大的权重上。"""
    if not weights:
        return weights
    delta = 1.0 - sum(weights)
    if delta:
        largest = max(range(len(weights)), key=lambda j: weights[j])
        weights[largest] += delta
    return weights


# ---------------------------------------------------------------------------
# 面板验证与提取
# ---------------------------------------------------------------------------

def validate(value):
    """验证平衡面板与配置，返回后续分析使用的规范化结构。"""
    if not isinstance(value, dict):
        raise ValueError("资料顶层必须是对象")

    weeks = value.get("weeks")
    if not isinstance(weeks, list) or not weeks:
        raise ValueError("weeks 必须是非空整数列表")
    for week in weeks:
        if not _is_int(week):
            raise ValueError("weeks 必须是连续整数")
    if weeks != list(range(weeks[0], weeks[-1] + 1)):
        raise ValueError("weeks 必须是连续整数周")

    sites = value.get("sites")
    if not isinstance(sites, list) or not sites:
        raise ValueError("sites 必须是非空厂区列表")
    seen_sites = set()
    for site in sites:
        if not isinstance(site, str) or not site or site in seen_sites:
            raise ValueError("厂区身份必须非空且唯一")
        seen_sites.add(site)

    rows = value.get("rows")
    if not isinstance(rows, list):
        raise ValueError("rows 必须是列表")
    week_set = set(weeks)
    values = {}
    week_index = {week: index for index, week in enumerate(weeks)}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("每行必须是对象")
        site = row.get("site")
        week = row.get("week")
        if not isinstance(site, str) or site not in seen_sites:
            raise ValueError("面板身份错误：未知厂区")
        if not _is_int(week) or week not in week_set:
            raise ValueError("面板身份错误：未知周")
        key = (site, week)
        if key in values:
            raise ValueError("面板身份错误：厂区/周重复")
        raw = row.get("value")
        if isinstance(raw, bool) or not isinstance(raw, (int, float)) or not math.isfinite(raw):
            raise ValueError("数值必须有限且不能是布尔值")
        values[key] = float(raw)
    if len(values) != len(sites) * len(weeks):
        raise ValueError("面板不是完整平衡资料")

    treated = value.get("treated")
    if not isinstance(treated, str) or treated not in seen_sites:
        raise ValueError("treated 必须是资料中的厂区")
    donors = [site for site in sites if site != treated]
    if not donors:
        raise ValueError("至少需要一个始终未处理的 donor")

    treatment_week = value.get("treatment_week")
    if not _is_int(treatment_week) or treatment_week not in week_set:
        raise ValueError("treatment_week 必须是面板中的周")
    split = week_index[treatment_week]
    pre_weeks = weeks[:split]
    post_weeks = weeks[split:]
    if not pre_weeks:
        raise ValueError("上线前至少需要一周")
    if not post_weeks:
        raise ValueError("上线后至少需要一周")

    validation_weeks = value.get("validation_weeks")
    if not _is_int(validation_weeks) or validation_weeks < 1:
        raise ValueError("validation_weeks 必须是正整数")
    if validation_weeks >= len(pre_weeks):
        raise ValueError("留出周必须非空且训练周必须非空")
    training_weeks = pre_weeks[:-validation_weeks]
    holdout_weeks = pre_weeks[-validation_weeks:]

    lambdas = value.get("lambdas")
    if not isinstance(lambdas, list) or not lambdas:
        raise ValueError("lambdas 必须是非空候选列表")
    lambda_values = []
    for index, lam in enumerate(lambdas):
        lambda_values.append(_number(lam, "lambdas[%d]" % index, nonnegative=True))

    pre_rmspe_max = _number(value.get("pre_rmspe_max"), "pre_rmspe_max", nonnegative=True)

    solver = value.get("solver", {})
    if not isinstance(solver, dict):
        raise ValueError("solver 必须是对象")
    tolerance = _number(solver.get("tolerance", DEFAULT_SOLVER_TOLERANCE), "solver.tolerance", positive=True)
    max_iterations = solver.get("max_iterations", DEFAULT_MAX_ITERATIONS)
    if not _is_int(max_iterations) or max_iterations < 1:
        raise ValueError("solver.max_iterations 必须是正整数")

    tie_tolerance = _number(value.get("tie_tolerance", DEFAULT_TIE_TOLERANCE), "tie_tolerance", nonnegative=True)
    zero_tolerance = _number(value.get("zero_tolerance", DEFAULT_ZERO_TOLERANCE), "zero_tolerance", nonnegative=True)

    series = {site: [values[(site, week)] for week in weeks] for site in sites}
    return {
        "weeks": weeks,
        "week_index": week_index,
        "sites": sites,
        "treated": treated,
        "donors": donors,
        "series": series,
        "treatment_week": treatment_week,
        "pre_weeks": pre_weeks,
        "post_weeks": post_weeks,
        "training_weeks": training_weeks,
        "holdout_weeks": holdout_weeks,
        "validation_weeks": validation_weeks,
        "lambdas": sorted(lambda_values),
        "pre_rmspe_max": pre_rmspe_max,
        "tie_tolerance": tie_tolerance,
        "zero_tolerance": zero_tolerance,
        "solver": {"tolerance": tolerance, "max_iterations": max_iterations},
    }


# ---------------------------------------------------------------------------
# 无矩阵求逆的单纯形约束二次规划
# ---------------------------------------------------------------------------

def _jacobi_eigen(matrix):
    """对称矩阵的 Jacobi 特征分解，返回 (特征值, 特征向量列)。

    只用初等旋转，不调用外部线性代数库；用来在共线/重复列造成奇异时
    求最小范数解，避免对 KKT 矩阵直接求逆。
    """
    n = len(matrix)
    a = [list(matrix[i]) for i in range(n)]
    # 消除输入中的微小不对称。
    for i in range(n):
        for j in range(i + 1, n):
            average = 0.5 * (a[i][j] + a[j][i])
            a[i][j] = a[j][i] = average
    vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    scale = max([abs(a[i][j]) for i in range(n) for j in range(n)] + [0.0])
    if scale == 0.0:
        return [0.0] * n, vectors
    for _ in range(100):
        off_diagonal = math.sqrt(sum(a[i][j] * a[i][j] for i in range(n) for j in range(i + 1, n)))
        if off_diagonal <= 1e-15 * scale:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) <= 1e-300:
                    continue
                angle = 0.5 * math.atan2(2.0 * apq, a[p][p] - a[q][q])
                cosine = math.cos(angle)
                sine = math.sin(angle)
                app = a[p][p]
                aqq = a[q][q]
                a[p][p] = cosine * cosine * app + 2.0 * cosine * sine * apq + sine * sine * aqq
                a[q][q] = sine * sine * app - 2.0 * cosine * sine * apq + cosine * cosine * aqq
                a[p][q] = a[q][p] = 0.0
                for k in range(n):
                    if k == p or k == q:
                        continue
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = a[p][k] = cosine * akp + sine * akq
                    a[k][q] = a[q][k] = -sine * akp + cosine * akq
                for k in range(n):
                    vkp = vectors[k][p]
                    vkq = vectors[k][q]
                    vectors[k][p] = cosine * vkp + sine * vkq
                    vectors[k][q] = -sine * vkp + cosine * vkq
    return [a[i][i] for i in range(n)], vectors


def _solve_linear(matrix, rhs):
    """用特征分解求对称半正定系统的最小范数解，奇异方向自动截断。"""
    n = len(rhs)
    if n == 0:
        return []
    eigenvalues, eigenvectors = _jacobi_eigen(matrix)
    largest = max(eigenvalues + [0.0])
    if largest <= 0.0:
        return [0.0] * n
    threshold = largest * 1e-13
    solution = [0.0] * n
    for index, eigenvalue in enumerate(eigenvalues):
        if eigenvalue <= threshold:
            continue
        projection = sum(eigenvectors[j][index] * rhs[j] for j in range(n))
        coefficient = projection / eigenvalue
        for j in range(n):
            solution[j] += coefficient * eigenvectors[j][index]
    return solution


def _face_candidate(hessian, linear, current, active):
    """在给定支持集上解等式约束 QP，返回候选向量（可能含负权重）。

    约束 sum(w)=1 通过 w = w0 + Z t 参数化；Z 的列张成 sum(z)=0。
    这样不需要对含约束的 KKT 矩阵直接求逆，共线方向也不会让流程失败。
    """
    n = len(current)
    k = len(active)
    candidate = [0.0] * n
    if k == 0:
        return candidate
    if k == 1:
        candidate[active[0]] = 1.0
        return candidate

    # w0 是支持集上的均匀可行点。
    w0 = [1.0 / k] * k
    # Z: k x (k-1) 的正交基，张成 sum(z)=0。使用 Helmert 基可以让
    # 奇异方向上的最小范数解等价于对称选择，重复列不会得到不对称权重。
    z = [[0.0] * (k - 1) for _ in range(k)]
    for col in range(k - 1):
        m = col + 1
        scale = 1.0 / math.sqrt(m * (m + 1))
        for row in range(m):
            z[row][col] = scale
        z[m][col] = -m * scale

    h_restricted = [[hessian[active[i]][active[j]] for j in range(k)] for i in range(k)]
    hw0 = [sum(h_restricted[i][j] * w0[j] for j in range(k)) for i in range(k)]
    c_restricted = [linear[active[i]] for i in range(k)]

    hz = [[sum(h_restricted[i][j] * z[j][col] for j in range(k)) for col in range(k - 1)] for i in range(k)]
    reduced = [
        [sum(z[i][row] * hz[i][col] for i in range(k)) for col in range(k - 1)]
        for row in range(k - 1)
    ]
    rhs = [
        -sum(z[i][col] * (hw0[i] + c_restricted[i]) for i in range(k))
        for col in range(k - 1)
    ]
    step = _solve_linear(reduced, rhs)
    for i, index in enumerate(active):
        candidate[index] = w0[i] + sum(z[i][col] * step[col] for col in range(k - 1))
    return candidate


def _kkt_report(hessian, linear, weights):
    """按项目约定的尺度给出 KKT 证据。

    L = f(w) + nu(sum(w)-1) - mu^T w，因此 mu = gradient + nu，
    mu >= 0 且 mu_j w_j = 0。等式乘子取支持集上梯度的平均值。
    """
    n = len(weights)
    gradient = [sum(hessian[i][j] * weights[j] for j in range(n)) + linear[i] for i in range(n)]
    support = [j for j, weight in enumerate(weights) if weight > ACTIVE_TOL]
    if support:
        multiplier = -sum(gradient[j] for j in support) / len(support)
    else:
        multiplier = -min(gradient) if gradient else 0.0
    mu = [gradient[j] + multiplier for j in range(n)]
    equality = abs(sum(weights) - 1.0)
    lower_violation = max(0.0, -min(weights)) if weights else 0.0
    active_stationarity = max([abs(mu[j]) for j in support] + [0.0])
    inactive = [j for j in range(n) if j not in support]
    dual_violation = max([max(0.0, -mu[j]) for j in inactive] + [0.0])
    complementarity = max([abs(mu[j] * weights[j]) for j in range(n)] + [0.0])
    maximum = max(equality, lower_violation, active_stationarity, dual_violation, complementarity)
    return {
        "equality_residual": equality,
        "lower_bound_violation": lower_violation,
        "active_stationarity_residual": active_stationarity,
        "dual_feasibility_violation": dual_violation,
        "nonnegativity_margin": min(mu) if mu else 0.0,
        "complementarity_residual": complementarity,
        "max_residual": maximum,
        "equality_multiplier": multiplier,
        "mu": mu,
        "gradient": gradient,
    }


def _solve_simplex_qp(hessian, linear, tolerance, max_iterations):
    """在 w>=0、sum(w)=1 上最小化 0.5*w'Hw + c'w。

    使用支持集/有效集方法：每次在支持集上解等式约束 QP；若候选越出
    非负边界，就沿可行方向走到边界并移除变量；若候选可行，则检查
    KKT，并把违反对偶可行性的变量加入支持集。不会显式求逆矩阵。
    """
    n = len(linear)
    weights = [1.0 / n] * n
    _fix_sum(weights)
    active = list(range(n))
    iterations = 0
    final_report = None

    for iteration in range(max_iterations):
        iterations = iteration + 1
        if not active:
            active = [max(range(n), key=lambda j: weights[j])]
        candidate = _face_candidate(hessian, linear, weights, active)
        minimum = min(candidate[j] for j in active)

        if minimum < -ACTIVE_TOL:
            # 向候选点走，直到第一个权重触到 0，然后缩小支持集。
            direction = [candidate[j] - weights[j] for j in range(n)]
            alpha = 1.0
            for j in active:
                if direction[j] < -1e-15:
                    alpha = min(alpha, weights[j] / (-direction[j]))
            alpha = max(0.0, min(1.0, alpha))
            weights = [weights[j] + alpha * direction[j] for j in range(n)]
            for j in range(n):
                if weights[j] < 0.0:
                    weights[j] = 0.0
            _fix_sum(weights)
            active = [j for j in range(n) if weights[j] > ACTIVE_TOL]
            if not active:
                active = [max(range(n), key=lambda j: weights[j])]
            continue

        # 候选可行：把它作为当前点，清掉数值上等于 0 的分量。
        weights = [max(0.0, candidate[j]) for j in range(n)]
        for j in range(n):
            if weights[j] <= ACTIVE_TOL:
                weights[j] = 0.0
        _fix_sum(weights)
        active = [j for j in range(n) if weights[j] > ACTIVE_TOL]
        if not active:
            active = [max(range(n), key=lambda j: weights[j])]
        final_report = _kkt_report(hessian, linear, weights)
        if final_report["max_residual"] <= tolerance:
            return {
                "weights": weights,
                "iterations": iterations,
                "converged": True,
                "termination_reason": "kkt_converged",
                "kkt": final_report,
            }

        inactive = [j for j in range(n) if j not in active]
        violating = [j for j in inactive if final_report["mu"][j] < 0.0]
        if not violating:
            # 理论上支持集内已最优；若仍超差，通常是数值尺度问题。
            return {
                "weights": weights,
                "iterations": iterations,
                "converged": False,
                "termination_reason": "stalled_without_dual_violation",
                "kkt": final_report,
            }
        join = min(violating, key=lambda j: final_report["mu"][j])
        active = sorted(active + [join])

    final_report = _kkt_report(hessian, linear, weights)
    return {
        "weights": weights,
        "iterations": iterations,
        "converged": False,
        "termination_reason": "budget_exhausted",
        "kkt": final_report,
    }


# ---------------------------------------------------------------------------
# 权重拟合与 lambda 选择
# ---------------------------------------------------------------------------

def fit_weights(donor_series, observed, lam, tolerance=DEFAULT_SOLVER_TOLERANCE,
                max_iterations=DEFAULT_MAX_ITERATIONS):
    """用 donor_series 拟合 observed，返回权重、目标与 KKT 证据。

    donor_series 是 J 条等长序列；observed 是同长度处理厂区序列。
    目标为 (1/T)||observed - synthetic||^2 + lam||w||^2，无截距。
    """
    lam = _number(lam, "lambda", nonnegative=True)
    tolerance = _number(tolerance, "solver.tolerance", positive=True)
    if not _is_int(max_iterations) or max_iterations < 1:
        raise ValueError("solver.max_iterations 必须是正整数")
    if not donor_series:
        raise ValueError("至少需要一个 donor")
    length = len(observed)
    if length == 0:
        raise ValueError("拟合周不能为空")
    for series in donor_series:
        if len(series) != length:
            raise ValueError("donor 与处理厂区序列长度不一致")
        for value in series:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError("数值必须有限且不能是布尔值")

    donors = [[float(value) for value in series] for series in donor_series]
    for value in observed:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("数值必须有限且不能是布尔值")
    target = [float(value) for value in observed]
    j = len(donors)
    # H = 2/T X'X + 2lambda I, c = -2/T X'y。
    hessian = [[0.0] * j for _ in range(j)]
    for row in range(j):
        for col in range(row, j):
            value = 2.0 / length * sum(donors[row][t] * donors[col][t] for t in range(length))
            if row == col:
                value += 2.0 * lam
            hessian[row][col] = value
            hessian[col][row] = value
    linear = [-2.0 / length * sum(donors[row][t] * target[t] for t in range(length)) for row in range(j)]

    solved = _solve_simplex_qp(hessian, linear, tolerance, max_iterations)
    weights = solved["weights"]
    synthetic = [sum(weights[j] * donors[j][t] for j in range(j)) for t in range(length)]
    residuals = [target[t] - synthetic[t] for t in range(length)]
    train_mse = sum(residual * residual for residual in residuals) / length
    objective = train_mse + lam * sum(weight * weight for weight in weights)

    kkt = solved["kkt"]
    return {
        "weights": weights,
        "gradient": kkt["gradient"],
        "equality_multiplier": kkt["equality_multiplier"],
        "mu": kkt["mu"],
        "kkt": {key: value for key, value in kkt.items() if key not in ("gradient", "mu", "equality_multiplier")},
        "objective": objective,
        "train_mse": train_mse,
        "train_rmspe": math.sqrt(train_mse),
        "iterations": solved["iterations"],
        "converged": solved["converged"],
        "optimal": bool(solved["converged"]),
        "termination_reason": solved["termination_reason"],
    }


def _fit_candidate(panel, donor_series, observed, lam):
    """在指定周集合上拟合一个单元，返回 fit 与留出预测。"""
    training = [panel["week_index"][week] for week in panel["training_weeks"]]
    holdout = [panel["week_index"][week] for week in panel["holdout_weeks"]]
    train_donors = [[series[index] for index in training] for series in donor_series]
    train_observed = [observed[index] for index in training]
    holdout_donors = [[series[index] for index in holdout] for series in donor_series]
    holdout_observed = [observed[index] for index in holdout]

    fit = fit_weights(
        train_donors,
        train_observed,
        lam,
        tolerance=panel["solver"]["tolerance"],
        max_iterations=panel["solver"]["max_iterations"],
    )
    holdout_mse = _mse(holdout_observed, fit["weights"], holdout_donors)
    return {
        "lambda": lam,
        "train_mse": fit["train_mse"],
        "train_rmspe": fit["train_rmspe"],
        "holdout_mse": holdout_mse,
        "objective": fit["objective"],
        "converged": fit["converged"],
        "optimal": fit["optimal"],
        "termination_reason": fit["termination_reason"],
        "iterations": fit["iterations"],
        "kkt_residual": fit["kkt"]["max_residual"],
        "kkt": fit["kkt"],
        "gradient": fit["gradient"],
        "equality_multiplier": fit["equality_multiplier"],
        "mu": fit["mu"],
        "weights": fit["weights"],
    }


def _select_lambda(panel, donor_series, observed):
    """按留出段均方误差选 lambda；同分（差值 <= tie_tolerance）取较大 lambda。"""
    candidates = [_fit_candidate(panel, donor_series, observed, lam) for lam in panel["lambdas"]]
    best = min(candidate["holdout_mse"] for candidate in candidates)
    tolerance = panel["tie_tolerance"]
    tied = [candidate for candidate in candidates if candidate["holdout_mse"] <= best + tolerance]
    selected = max(tied, key=lambda candidate: candidate["lambda"])
    return {
        "selected": selected["lambda"],
        "best_holdout_mse": best,
        "tie_tolerance": tolerance,
        "selected_holdout_mse": selected["holdout_mse"],
        "selection_rule": "min holdout MSE; absolute tie within tie_tolerance -> larger lambda",
        "candidates": candidates,
    }


def _fit_protocol(panel, donor_series, observed):
    """完整协议：选 lambda（仅训练与留出段），再在全部上线前周重拟合。"""
    selection = _select_lambda(panel, donor_series, observed)
    pre_indices = [panel["week_index"][week] for week in panel["pre_weeks"]]
    pre_donors = [[series[index] for index in pre_indices] for series in donor_series]
    pre_observed = [observed[index] for index in pre_indices]
    final = fit_weights(
        pre_donors,
        pre_observed,
        selection["selected"],
        tolerance=panel["solver"]["tolerance"],
        max_iterations=panel["solver"]["max_iterations"],
    )
    return selection, final


# ---------------------------------------------------------------------------
# 比值、差距与诊断输出
# ---------------------------------------------------------------------------

def _gap_details(panel, weights, observed, donor_series):
    weekly = []
    for index, week in enumerate(panel["weeks"]):
        synthetic = sum(weight * series[index] for weight, series in zip(weights, donor_series))
        weekly.append({
            "week": week,
            "period": "pre" if week < panel["treatment_week"] else "post",
            "observed": observed[index],
            "synthetic": synthetic,
            "gap": observed[index] - synthetic,
        })
    pre_gaps = [entry["gap"] for entry in weekly if entry["period"] == "pre"]
    post_gaps = [entry["gap"] for entry in weekly if entry["period"] == "post"]
    pre_msp = sum(gap * gap for gap in pre_gaps) / len(pre_gaps)
    post_msp = sum(gap * gap for gap in post_gaps) / len(post_gaps)
    return {
        "weekly": weekly,
        "pre_msp": pre_msp,
        "pre_rmspe": math.sqrt(pre_msp),
        "post_msp": post_msp,
        "post_rmse": math.sqrt(post_msp),
        "post_average_gap": sum(post_gaps) / len(post_gaps),
    }


def _ratio_report(pre_msp, post_msp, zero_tolerance):
    """按项目约定处理 0/0 与 正/0，绝不输出 NaN/Infinity 常量。"""
    if pre_msp <= zero_tolerance:
        if post_msp <= zero_tolerance:
            return {
                "same_period_ratio": 0.0,
                "ratio_kind": "zero_over_zero",
                "ratio_display": 0.0,
                "pre_msp": pre_msp,
                "post_msp": post_msp,
            }
        return {
            "same_period_ratio": float("inf"),
            "ratio_kind": "positive_over_zero",
            "ratio_display": _POSITIVE_INFINITY,
            "pre_msp": pre_msp,
            "post_msp": post_msp,
        }
    return {
        "same_period_ratio": post_msp / pre_msp,
        "ratio_kind": "finite",
        "ratio_display": post_msp / pre_msp,
        "pre_msp": pre_msp,
        "post_msp": post_msp,
    }


def _ratio_equal(left, right, tolerance):
    if left == right:
        return True
    if math.isinf(left) or math.isinf(right):
        return False
    return abs(left - right) <= tolerance


# ---------------------------------------------------------------------------
# 主分析与空间 placebo 诊断
# ---------------------------------------------------------------------------

def analyze(value, source_sha256=None):
    """在平衡周面板上执行主合成对照与空间 placebo 诊断。

    value 可以是已解析的 JSON 对象，也可以是资料路径；路径模式会自动
    计算原始文件指纹。
    """
    if isinstance(value, (str, Path)):
        value, source_sha256 = load(value)
    panel = validate(value)
    donors = panel["donors"]
    donor_series = [panel["series"][donor] for donor in donors]
    observed = panel["series"][panel["treated"]]

    selection, final = _fit_protocol(panel, donor_series, observed)
    details = _gap_details(panel, final["weights"], observed, donor_series)
    main_ratio = _ratio_report(details["pre_msp"], details["post_msp"], panel["zero_tolerance"])

    pre_fit_ok = details["pre_rmspe"] <= panel["pre_rmspe_max"]
    warning = None
    if not pre_fit_ok:
        warning = (
            "上线前拟合失配：pre-RMSPE=%.12g 超过门槛 %.12g。结果仍照常展示，"
            "但该处理厂区在配置的 placebo 门槛下不可比；不得用处理后结果挑选 donor 或 lambda。"
            % (details["pre_rmspe"], panel["pre_rmspe_max"])
        )

    result = {
        "source_sha256": source_sha256,
        "protocol": {
            "treated": panel["treated"],
            "donors": donors,
            "weeks": panel["weeks"],
            "treatment_week": panel["treatment_week"],
            "pre_weeks": panel["pre_weeks"],
            "post_weeks": panel["post_weeks"],
            "training_weeks": panel["training_weeks"],
            "holdout_weeks": panel["holdout_weeks"],
            "validation_weeks": panel["validation_weeks"],
            "lambdas": panel["lambdas"],
            "tie_tolerance": panel["tie_tolerance"],
            "zero_tolerance": panel["zero_tolerance"],
            "pre_rmspe_max": panel["pre_rmspe_max"],
            "solver": panel["solver"],
            "intercept": False,
            "standardization": "none",
            "constraints": "w>=0, sum(w)=1",
            "objective": "training MSE + lambda * sum(w^2)",
        },
        "lambda_selection": {
            "selected": selection["selected"],
            "best_holdout_mse": selection["best_holdout_mse"],
            "selected_holdout_mse": selection["selected_holdout_mse"],
            "tie_tolerance": selection["tie_tolerance"],
            "selection_rule": selection["selection_rule"],
            "candidates": selection["candidates"],
        },
        "fit": {
            "lambda": selection["selected"],
            "fit_weeks": panel["pre_weeks"],
            "fit_scope": "all pre weeks (refit after lambda selection)",
            "weights": [
                {"donor": donor, "weight": weight}
                for donor, weight in zip(donors, final["weights"])
            ],
            "donor_weights": {
                donor: weight for donor, weight in zip(donors, final["weights"])
            },
            "weight_sum": sum(final["weights"]),
            "objective": final["objective"],
            "training_mse": final["train_mse"],
            "training_rmspe": final["train_rmspe"],
            "gradient": final["gradient"],
            "equality_multiplier": final["equality_multiplier"],
            "mu": final["mu"],
            "pre_msp": details["pre_msp"],
            "pre_rmspe": details["pre_rmspe"],
            "pre_msp_is_zero": details["pre_msp"] <= panel["zero_tolerance"],
            "post_msp": details["post_msp"],
            "post_rmse": details["post_rmse"],
            "post_msp_is_zero": details["post_msp"] <= panel["zero_tolerance"],
            "post_average_gap": details["post_average_gap"],
            "average_post_gap": details["post_average_gap"],
            "same_period_ratio": main_ratio["ratio_display"],
            "post_msp_over_pre_msp": main_ratio["ratio_display"],
            "ratio": main_ratio["ratio_display"],
            "ratio_kind": main_ratio["ratio_kind"],
            "main_statistic": {
                "name": "post_msp_over_pre_msp",
                "pre_msp": details["pre_msp"],
                "post_msp": details["post_msp"],
                "ratio": main_ratio["ratio_display"],
                "ratio_kind": main_ratio["ratio_kind"],
            },
            "kkt": final["kkt"],
            "solver": {
                "converged": final["converged"],
                "optimal": final["optimal"],
                "termination_reason": final["termination_reason"],
                "iterations": final["iterations"],
            },
        },
        "weekly": details["weekly"],
        "post_gap_curve": [
            entry for entry in details["weekly"] if entry["period"] == "post"
        ],
        "diagnostics": {
            "pre_fit_within_threshold": pre_fit_ok,
            "mismatch_warning": warning,
            "causal_disclaimer": (
                "合成对照按可观测的 donor 周序列构造反事实。它不能自动消除未观察混杂；"
                "差距下降是拟合权重下的对照结果，不是投入有效的因果证明。"
            ),
        },
    }
    result["spatial_placebo"] = _spatial_placebo(
        panel,
        main_ratio["same_period_ratio"],
        main_ratio["ratio_display"],
        main_ratio["ratio_kind"],
        details["pre_msp"],
        details["pre_rmspe"],
        pre_fit_ok,
    )
    return result


def _placebo_unit(panel, pseudo, donor_pool):
    """按同一协议把一个 donor 当作伪处理厂区。"""
    observed = panel["series"][pseudo]
    donor_series = [panel["series"][donor] for donor in donor_pool]
    selection, final = _fit_protocol(panel, donor_series, observed)
    details = _gap_details(panel, final["weights"], observed, donor_series)
    ratio = _ratio_report(details["pre_msp"], details["post_msp"], panel["zero_tolerance"])
    eligible = details["pre_rmspe"] <= panel["pre_rmspe_max"]
    exclusion_reason = None if eligible else "pre_rmspe_exceeds_threshold"
    return {
        "donor": pseudo,
        "donor_pool": donor_pool,
        "lambda_selected": selection["selected"],
        "holdout_mse": selection["selected_holdout_mse"],
        "pre_msp": details["pre_msp"],
        "pre_rmspe": details["pre_rmspe"],
        "pre_msp_is_zero": details["pre_msp"] <= panel["zero_tolerance"],
        "post_msp": details["post_msp"],
        "post_msp_is_zero": details["post_msp"] <= panel["zero_tolerance"],
        "post_average_gap": details["post_average_gap"],
        "same_period_ratio": ratio["ratio_display"],
        "post_msp_over_pre_msp": ratio["ratio_display"],
        "ratio": ratio["ratio_display"],
        "ratio_kind": ratio["ratio_kind"],
        "_ratio_value": ratio["same_period_ratio"],
        "pre_fit_within_threshold": eligible,
        "eligible": eligible,
        "exclusion_reason": exclusion_reason,
        "weights": [
            {"donor": donor, "weight": weight}
            for donor, weight in zip(donor_pool, final["weights"])
        ],
        "solver": {
            "converged": final["converged"],
            "optimal": final["optimal"],
            "termination_reason": final["termination_reason"],
            "iterations": final["iterations"],
            "kkt_residual": final["kkt"]["max_residual"],
        },
    }


def _spatial_placebo(panel, main_ratio_value, main_ratio_display, main_ratio_kind,
                     main_pre_msp, main_pre_rmspe, main_pre_ok):
    """空间 placebo 诊断：轮换伪处理厂区，原受处理厂区不进任何 donor 池。"""
    placebos = []
    for donor in panel["donors"]:
        pool = [other for other in panel["donors"] if other != donor]
        if not pool:
            placebos.append({
                "donor": donor,
                "donor_pool": [],
                "eligible": False,
                "exclusion_reason": "empty_placebo_donor_pool",
                "same_period_ratio": None,
                "ratio_kind": None,
                "pre_rmspe": None,
                "pre_msp": None,
                "post_msp": None,
                "_ratio_value": None,
                "weights": [],
                "solver": None,
            })
            continue
        placebos.append(_placebo_unit(panel, donor, pool))

    eligible = [entry for entry in placebos if entry["eligible"]]
    excluded = [
        {
            "donor": entry["donor"],
            "pre_rmspe": entry["pre_rmspe"],
            "threshold": panel["pre_rmspe_max"],
            "reason": entry["exclusion_reason"],
        }
        for entry in placebos if not entry["eligible"]
    ]

    if not main_pre_ok:
        rank = {
            "reportable": False,
            "reason": "treated_pre_rmspe_exceeds_threshold",
            "treated_pre_rmspe": main_pre_rmspe,
            "treated_ratio": main_ratio_display,
            "treated_ratio_kind": main_ratio_kind,
            "threshold": panel["pre_rmspe_max"],
            "denominator": None,
            "numerator": None,
            "proportion": None,
            "tail_proportion": None,
        }
        tied = []
    else:
        denominator = 1 + len(eligible)
        numerator = 1 + sum(
            1 for entry in eligible
            if entry["_ratio_value"] is not None and entry["_ratio_value"] >= main_ratio_value
        )
        tied = [
            entry["donor"] for entry in eligible
            if entry["_ratio_value"] is not None
            and _ratio_equal(entry["_ratio_value"], main_ratio_value, panel["tie_tolerance"])
        ]
        rank = {
            "reportable": True,
            "reason": None,
            "treated_pre_rmspe": main_pre_rmspe,
            "treated_ratio": main_ratio_display,
            "treated_ratio_kind": main_ratio_kind,
            "threshold": panel["pre_rmspe_max"],
            "denominator": denominator,
            "numerator": numerator,
            "proportion": numerator / denominator,
            "tail_proportion": numerator / denominator,
            "tie_rule": "ratio >= treated ratio; ties included",
        }

    # 去掉仅用于比较的内部比值。
    public_placebos = []
    for entry in placebos:
        entry = dict(entry)
        entry.pop("_ratio_value", None)
        public_placebos.append(entry)

    return {
        "name": "空间 placebo 诊断",
        "not_a_randomization_p_value": True,
        "original_treated_excluded_from_placebo_pools": True,
        "placebo_pool_rule": "for each donor, pool = all donors except that donor; original treated is never a donor",
        "disclaimer": (
            "这是空间 placebo 诊断，不是真实随机分配下的精确 p 值；"
            "它只描述受处理厂区在本组 donor 轮换中的相对位置。"
        ),
        "threshold": panel["pre_rmspe_max"],
        "treated": {
            "pre_msp": main_pre_msp,
            "pre_rmspe": main_pre_rmspe,
            "pre_fit_within_threshold": main_pre_ok,
            "same_period_ratio": main_ratio_display,
            "ratio_kind": main_ratio_kind,
        },
        "eligible_placebo_count": len(eligible),
        "included_denominator": rank["denominator"],
        "included_numerator": rank["numerator"],
        "tail_rank_proportion": rank["proportion"],
        "rank": rank,
        "ties_included_in_numerator": True,
        "tied_placebos": tied,
        "excluded": excluded,
        "placebos": public_placebos,
    }
