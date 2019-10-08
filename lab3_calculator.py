# Griffin Saiia, gjs64
# E&M Lab 3: EOM
# Data processing
import os
import sys
from math import sqrt, pi

data_folder="data/"
error = False

# constants
u=(4*pi)*(10**-7)
n=130
r_coil=0.158
uncert_r = 0.005

uncert_distance = 0.05*(10**-2)
uncert_powerSupply = 0.001

# for constant voltage
i_slope = 0.0
i_intercept = 0.0

# for constant current
v_slope = 0.0
v_intercept = 0.0

# calculated e/m
v_const = 0.0
i_e_div_m = 0.0
best_i = ["", "", "", ""]
i_slope_e_div_m = (0.0, 0.0)

i_const = 0.0
v_e_div_m = 0.0
best_v = ["", "", "", ""]
v_slope_e_div_m = (0.0, 0.0)

# expected e/m
e_div_m=1.75882017*(10**11)
e_div_m_unc=0.00000007*(10**11)

def process_iTuples(const, i_vs_R):
    global i_e_div_m, i_slope_e_div_m, best_i, v_const, i_slope, i_intercept, i_slope_uncert
    v_const = const[0]
    i_slope = const[1]
    i_slope_uncert = const[3]
    i_intercept = const[2]
    em_Values = []
    em_Uncert = []
    for each in i_vs_R:
        i = each[0]
        r = each[1]
        b = calcMagField(i)
        b_uncert = uncertMagField(i)
        em_Values.append(calc_e_div_m(v_const,b,r))
        em_Uncert.append(uncert_e_div_m(v_const, b, b_uncert, r))
    # average e/m value according to data
    i_e_div_m = ((float(sum(em_Values))/len(em_Values)), (float(sum(em_Uncert)/len(em_Uncert))))
    i_slope_e_div_m = (calcEm_wISlope(i_slope, v_const), uncertEm_wISlope(i_slope, i_slope_uncert, v_const, (v_const * uncert_powerSupply)))
    min = getBestEmValue(em_Values)
    b_best = calcMagField(i_vs_R[min[0]][0])
    b_uncert_best = uncertMagField(i_vs_R[min[0]][0])
    em_uncert_best = uncert_e_div_m(v_const, b_best, b_uncert_best, i_vs_R[min[0]][1])
    best_i = [str(i_vs_R[min[0]][0])+" A", str(min[1]), str(checkSig(min[2], em_uncert_best)), str(em_uncert_best)]

def process_vTuples(const, v_vs_R):
    global v_e_div_m, v_slope_e_div_m, best_v, i_const, v_slope, v_intercept, v_slope_uncert
    i_const = const[0]
    v_slope = const[1]
    v_slope_uncert = const[3]
    v_intercept = const[2]
    em_Values = []
    em_Uncert = []
    for each in v_vs_R:
        v = each[0]
        r = each[1]
        b = calcMagField(i_const)
        b_uncert = uncertMagField(i_const)
        em_Values.append(calc_e_div_m(v,b,r))
        em_Uncert.append(uncert_e_div_m(v,b,b_uncert,r))
    # average e/m value according to data
    v_e_div_m = ((float(sum(em_Values))/len(em_Values)), (float(sum(em_Uncert))/len(em_Uncert)))
    v_slope_e_div_m = (calcEm_wVSlope(v_slope, i_const), uncertEm_wVSlope(v_slope, v_slope_uncert, i_const, (i_const*uncert_powerSupply)))
    min = getBestEmValue(em_Values)
    b_best = calcMagField(i_const)
    b_uncert_best = uncertMagField(i_const)
    em_uncert_best = uncert_e_div_m(v_vs_R[min[0]][0], b_best, b_uncert_best, v_vs_R[min[0]][1])
    best_v = [str(v_vs_R[min[0]][0])+" V", str(min[1]), str(checkSig(min[2], em_uncert_best)), str(em_uncert_best)]

def printError(const, tuples):
    error = True
    print("Error stripping data:")
    print(" Consts (const, slope, intercept):")
    for each in const:
        print("     "+str(each))
    print(" DATA:")
    for each in tuples:
        print("     "+str(each))

func_map = {"V":process_iTuples,
            "I":process_vTuples,
            "ERROR":printError}

def main():
    files = os.listdir(data_folder)
    for each in files:
        (unit, const, tuples)=process_data(data_folder+"/"+each)
        func_map[unit](const, tuples)
    if not error:
        analyzeResult()

def process_data(file):
    lines = []
    const = []
    unit = "ERROR"
    tuples = []
    with open(file) as f:
        lines = f.readlines()
        f.close()
    i = 0
    try:
        for line in lines:
            if '#' not in line:
                split = line.split(",")
                if(i == 0):
                    const.append(float(split[1]))
                    unit = split[0]
                elif(i < 4):
                    const.append(float(split[1]))
                else:
                    meters = float(split[1])*(10**-2)
                    tuples.append( (float(split[0]), meters) )
                i += 1
        return (unit, const, tuples)
    except IndexError:
        tuples.append(("UNIT", unit))
        unit = "ERROR"
        return (unit, const, tuples)

def calcEm_wISlope(a, v):
    return (125*(a**2)*(r_coil**2)*(v))/(32*(n**2)*(u**2))

def uncertEm_wISlope(a, uncert_a, v, uncert_v):
    uncert_EM_r = ((250*(a**2)*r_coil*v)/(32*(n**2)*(u**2)))*uncert_r
    uncert_EM_a = ((250*a*(r_coil**2)*v)/(32*(n**2)*(u**2)))*uncert_a
    uncert_EM_v = ((125*(a**2)*(r_coil**2))/(32*(n**2)*(u**2)))*uncert_v
    return sqrt((uncert_EM_r**2)+(uncert_EM_a**2)+(uncert_EM_v**2))

def calcEm_wVSlope(b, i):
    return (125*(r_coil**2))/(32*(b**2)*(i**2)*(n**2)*(u**2))

def uncertEm_wVSlope(b, uncert_b, i, uncert_i):
    uncert_EM_r = ((250*r_coil)/(32*(b**2)*(i**2)*(n**2)*(u**2)))*uncert_r
    uncert_EM_b = ((250*(r_coil**2))/(32*(b**3)*(i**2)*(n**2)*(u**2)))*uncert_b
    uncert_EM_i = ((250*(r_coil**2))/(32*(b**2)*(i**3)*(n**2)*(u**2)))*uncert_i
    return sqrt((uncert_EM_r**2)+(uncert_EM_b**2)+(uncert_EM_i**2))

def calcMagField(i):
    return (8*u*n*i)/(5*r_coil*sqrt(5))

def uncertMagField(i):
    uncert_Mag_i = ((8*u*n)/(5*r_coil*sqrt(5)))*(i*uncert_powerSupply)
    uncert_Mag_r = ((8*u*n*i)/(5*(r_coil**2)*sqrt(5)))*uncert_r
    return sqrt((uncert_Mag_i**2)+(uncert_Mag_r**2))

def calc_e_div_m(v,b,r):
    return (2*v)/((b*r)**2)

def uncert_e_div_m(v,b,uncert_b,r):
    uncert_Em_v = (2/((b*r)**2))*(v*uncert_powerSupply)
    uncert_Em_b = ((4*v)/((b**3)*(r**2)))*uncert_b
    uncert_Em_r = ((4*v)/((b**2)*(r**3)))*uncert_r
    return sqrt((uncert_Em_v**2)+(uncert_Em_b**2)+(uncert_Em_r**2))

def diffExpectedCalculated(calc):
    if(calc > e_div_m):
        return (calc - e_div_m)
    else:
        return (e_div_m - calc)

def getBestEmValue(em_Values):
    i = 0
    diff = diffExpectedCalculated(em_Values[i])
    min = [i, em_Values[i], diff]
    for value in em_Values:
        diff = diffExpectedCalculated(value)
        if(diff < min[2]):
            min = [i, value, diff]
        i+=1
    return min

def checkSig(diff, uncert):
    return (diff < (uncert+e_div_m_unc))

def analyzeResult():
    print("With Varied Current, and Constant Voltage @ "+str(v_const)+" V:")
    printGraphCalc(i_slope, i_slope_uncert, i_intercept, i_slope_e_div_m[0], i_slope_e_div_m[1])
    printProgCalc(i_e_div_m[0], i_e_div_m[1], best_i)
    print("\n")
    print("With Varied Voltage, and Constant Current @ "+str(i_const)+" A:")
    printGraphCalc(v_slope, v_slope_uncert, v_intercept, v_slope_e_div_m[0], v_slope_e_div_m[1])
    printProgCalc(v_e_div_m[0], v_e_div_m[1], best_v)

def printGraphCalc(s, s_uncert, int, em, em_uncert):
    print(" Graphical Calculations:")
    print("     slope: "+str(s)+" +/- "+str(s_uncert))
    print("     intercept: "+str(int))
    print("     slope-based e/m calculation: "+str(em)+" +/- "+str(em_uncert))
    diff = diffExpectedCalculated(em)
    sig = checkSig(diff, em_uncert)
    print("     expected e/m value within slope-based error: "+str(sig))

def printProgCalc(em, em_uncert, best_val):
    print(" Program Analysis:")
    print("     average e/m value from data: "+str(em)+" +/- "+str(em_uncert))
    splitt = best_val[0].split(" ")
    justNum = splitt[0]
    unit = splitt[1]
    print("     best input for calculated e/m value: "+justNum+" +/- "+str(float(justNum) * uncert_powerSupply)+" "+unit)
    print("     e/m value for best input: "+best_val[1]+" +/- "+best_val[3])
    print("     expected e/m value within best e/m value error: "+best_val[2])

# Execute the wrapper
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print
        print 'Interrupted \_[o_0]_/'
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
