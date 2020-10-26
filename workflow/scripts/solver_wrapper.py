import pulp
import warnings
import sys

mip_solvs = ["SCIP_CMD", "CHOCO_CMD", "COINMP_DLL", "CPLEX_DLL", "GLPK_CMD", "LpSolver", "MIPCL_CMD",
             "MOSEK", "PULP_CHOCO_CMD", "PYGLPK", "YAPOSIB", "LpSolver_CMD", "CPLEX_PY", "GUROBI", "XPRESS",
             "GUROBI_CMD", "PULP_CBC_CMD", "COIN_CMD", "CPLEX_CMD"]
timelimit_solvs = mip_solvs[1:]
gapRel_solvs = mip_solvs[12:]
thread_solvs = mip_solvs[15:]


def create_solver(solver, mip, timelimit, gaprel, gapabs, fracgap, maxnodes, maxmemory, threads):
    if solver == "SCIP_CMD":
        return pulp.apis.SCIP_CMD(mip=mip)
    if solver == "CHOCO_CMD":
        return pulp.apis.CHOCO_CMD(mip=mip, timeLimit=timelimit)
    if solver == "COINMP_DLL":
        return pulp.apis.COINMP_DLL(mip=mip, timeLimit=timelimit)
    if solver == "CPLEX_DLL":
        return pulp.apis.CPLEX_DLL(mip=mip, timeLimit=timelimit)
    if solver == "GLPK_CMD":
        return pulp.apis.GLPK_CMD(mip=mip, timeLimit=timelimit)
    if solver == "LpSolver":
        return pulp.apis.LpSolver(mip=mip, timeLimit=timelimit)
    if solver == "MIPCL_CMD":
        return pulp.apis.MIPCL_CMD(mip=mip, timeLimit=timelimit)
    if solver == "MOSEK":
        return pulp.apis.MOSEK(mip=mip, timeLimit=timelimit)
    if solver == "PULP_CHOCO_CMD":
        return pulp.apis.PULP_CHOCO_CMD(mip=mip, timeLimit=timelimit)
    if solver == "PYGLPK":
        return pulp.apis.PYGLPK(mip=mip, timeLimit=timelimit)
    if solver == "YAPOSIB":
        return pulp.apis.YAPOSIB(mip=mip, timeLimit=timelimit)
    if solver == "LpSolver_CMD":
        return pulp.apis.LpSolver_CMD(mip=mip, timeLimit=timelimit)
    if solver == "CPLEX_PY":
        return pulp.apis.CPLEX_PY(mip=mip, timeLimit=timelimit, gapRel=gaprel)
    if solver == "GUROBI":
        return pulp.apis.GUROBI(mip=mip, timeLimit=timelimit, gapRel=gaprel)
    if solver == "XPRESS":
        return pulp.apis.XPRESS(mip=mip, timeLimit=timelimit, gapRel=gaprel)
    if solver == "GUROBI_CMD":
        return pulp.apis.GUROBI_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs)
    if solver == "PULP_CBC_CMD":
        return pulp.apis.PULP_CBC_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs)
    if solver == "COIN_CMD":
        return pulp.apis.COIN_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs,
                                  fracGap=fracgap)
    if solver == "CPLEX_CMD":
        return pulp.apis.CPLEX_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs,
                                   maxMemory=maxmemory, maxNodes=maxnodes)
    return ""


def check_params(solver, mip, timelimit, gaprel, gapabs, fracgap, maxnodes, maxmemory, threads):
    if not threads or threads < 1:
        threads = None
    if not mip:
        mip = False
    if mip:
        if not isinstance(mip, bool):
            warnings.warn("The mip option can only be set to True or False. The option was set to False.")
            mip = False
        if not solver in mip_solvs:
            warnings.warn("The mip option is not available for {}.".format(solver))
    if timelimit:
        if not isinstance(timelimit, int) or timelimit < 1:
            warnings.warn("The timeLimit option must be a positive integer number. The option was set to None.")
            timelimit = None
        if not solver in timelimit_solvs:
            warnings.warn("The timeLimit option is not available for {}.".format(solver))
    if gaprel:
        if not isinstance(gaprel, float):
            warnings.warn("The gapRel option must be a floating point number. The option was set to None.")
            gaprel = None
        if not solver in gapRel_solvs:
            warnings.warn("The gapRel option is not available for {}.".format(solver))
    if threads or gapabs:
        if not isinstance(threads, int):
            warnings.warn("The thread option must be a positive integer number. The option was set to None.")
            threads = None
        if not isinstance(gapabs, float):
            warnings.warn("The gapAbs option must be a floating point number. The option was set to None.")
            gapabs = None
        if not solver in thread_solvs:
            warnings.warn("The gapAbs or threads options are not available for {}.".format(solver))
    if fracgap:
        if not isinstance(fracgap, float):
            warnings.warn("The fracGap option must be a floating point number. The option was set to None.")
            fracgap = None
        if not solver == "COIN_CMD":
            warnings.warn("The fracGap option is not available for {}.".format(solver))
    if maxnodes or maxmemory:
        if not isinstance(maxnodes, int) or maxnodes < 1:
            warnings.warn("The maxNodes option must be a positive integer number. The option was set to None.")
            maxnodes = None
        if not isinstance(maxmemory, int) or maxmemory < 1:
            warnings.warn("The maxMemory option must be a positive integer number. The option was set to None.")
            maxmemory = None
        if not solver == "CPLEX_CMD":
            warnings.warn("The fracGap option is not available for {}.".format(solver))
    create_solver(solver, mip, timelimit, gaprel, gapabs, fracgap, maxnodes, maxmemory, threads)
