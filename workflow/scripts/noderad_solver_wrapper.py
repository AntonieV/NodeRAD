import pulp
import sys

mip_solvs = ["SCIP_CMD", "CHOCO_CMD", "COINMP_DLL", "CPLEX_DLL", "GLPK_CMD", "LpSolver", "MIPCL_CMD",
             "MOSEK", "PULP_CHOCO_CMD", "PYGLPK", "YAPOSIB", "LpSolver_CMD", "CPLEX_PY", "GUROBI", "XPRESS",
             "GUROBI_CMD", "PULP_CBC_CMD", "COIN_CMD", "CPLEX_CMD"]
timelimit_solvs = mip_solvs[1:]
gapRel_solvs = mip_solvs[12:]
thread_solvs = mip_solvs[15:]


def create_solver(solver, mip, timelimit, gaprel, gapabs, maxnodes, maxmemory, threads):
    sys.stderr.write("Params for ILP-solver: solver={} mip={} timeLimit={} gapRel={} gapAbs={} maxNodes={} maxMemory={} threads={}\n".format(str(solver), str(mip), str(timelimit), str(gaprel), str(gapabs), str(maxnodes), str(maxmemory), str(threads)))
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
        return pulp.apis.COIN_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs)
    if solver == "CPLEX_CMD":
        return pulp.apis.CPLEX_CMD(mip=mip, timeLimit=timelimit, threads=threads, gapRel=gaprel, gapAbs=gapabs,
                                   maxMemory=maxmemory, maxNodes=maxnodes)
    return ""


def check_params(solver, mip, timelimit, gaprel, gapabs, maxnodes, maxmemory, threads):
    if not threads or threads < 1:
        threads = None
    if not mip:
        mip = False
    if mip:
        if not isinstance(mip, bool):
            sys.stderr.write("The mip option can only be set to True or False, was setting to False.\n")
            mip = False
        if not solver in mip_solvs:
            sys.stderr.write("The mip option is not available for {}\n.".format(solver))
    if timelimit:
        if not isinstance(timelimit, int) or timelimit < 1:
            sys.stderr.write("The timeLimit option must be a positive integer number, setting to None.\n")
            timelimit = None
        if not solver in timelimit_solvs:
            sys.stderr.write("The timeLimit option is not available for {}.\n".format(solver))
    if gaprel:
        if not isinstance(gaprel, float):
            sys.stderr.write("The gapRel option must be a floating point number, setting to None.\n")
            gaprel = None
        if not solver in gapRel_solvs:
            sys.stderr.write("The gapRel option is not available for {}.\n".format(solver))
    if threads or gapabs:
        if not isinstance(threads, int):
            sys.stderr.write("The thread option must be a positive integer number, setting to None.\n")
            threads = None
        if not isinstance(gapabs, float):
            sys.stderr.write("The gapAbs option must be a floating point number, setting to None.\n")
            gapabs = None
        if not solver in thread_solvs:
            sys.stderr.write("The gapAbs or threads options are not available for {}.\n".format(solver))
    if maxnodes or maxmemory:
        if not isinstance(maxnodes, int) or maxnodes < 1:
            sys.stderr.write("The maxNodes option must be a positive integer number, setting to None.\n")
            maxnodes = None
        if not isinstance(maxmemory, int) or maxmemory < 1:
            sys.stderr.write("The maxMemory option must be a positive integer number, setting to None.\n")
            maxmemory = None
        if not solver == "CPLEX_CMD":
            sys.stderr.write("The fracGap option is not available for {}\n.".format(solver))
    create_solver(solver, mip, timelimit, gaprel, gapabs, maxnodes, maxmemory, threads)
