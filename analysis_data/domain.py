"""核对平衡周面板并查看各厂区原始均值。"""
from .io import finite
def summarize(value):
    weeks=value['weeks'];sites=value['sites'];rows=value['rows'];seen={}
    if not weeks or weeks!=list(range(weeks[0],weeks[-1]+1)) or len(set(sites))!=len(sites):raise ValueError('时间轴或厂区重复')
    for row in rows:
        key=(row['site'],row['week'])
        if key in seen or key[0] not in sites or key[1] not in weeks:raise ValueError('面板身份错误')
        seen[key]=finite(row['value'])
    if len(seen)!=len(sites)*len(weeks):raise ValueError('面板不是完整平衡资料')
    return [{'site':s,'mean':sum(seen[s,t] for t in weeks)/len(weeks),'weeks':len(weeks)} for s in sorted(sites)]
