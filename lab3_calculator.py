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

# for constant voltage
i_slope = 0.0
i_intercept = 0.0

# for constant current
v_slope = 0.0
v_intercept = 0.0

# calculated e/m
v_const = 0.0
i_e_div_m = 0.0
best_i = ["", 0.0, False]
i_slope_e_div_m = 0.0

i_const = 0.0
v_e_div_m = 0.0
best_v = ["", 0.0, False]
v_slope_e_div_m = 0.0

# expected e/m
e_div_m=1.75882017*(10**11)
e_div_m_unc=0.00000007*(10**11)

def process_iTuples(const, i_vs_R):
    global i_e_div_m, i_slope_e_div_m, best_i, v_const, i_slope, i_intercept
    v_const = const[0]
    i_slope = const[1]
    i_intercept = const[2]
    em_Values = []
    for each in i_vs_R:
        i = each[0]
        r = each[1]
        b = calcMagField(i)
        em_Values.append(calc_e_div_m(v_const,b,r))
    # average e/m value according to data
    i_e_div_m = float(sum(em_Values))/len(em_Values)
    i_slope_e_div_m = calcEm_wISlope(i_slope, v_const)
    min = getBestEmValue(em_Values)
    best_i = [str(i_vs_R[min[0]][0])+" A", str(min[1]), str(min[3])]

def process_vTuples(const, v_vs_R):
    global v_e_div_m, v_slope_e_div_m, best_v, i_const, v_slope, v_intercept
    i_const = const[0]
    v_slope = const[1]
    v_intercept = const[2]
    em_Values = []
    for each in v_vs_R:
        v = each[0]
        r = each[1]
        b = calcMagField(i_const)
        em_Values.append(calc_e_div_m(v,b,r))
    # average e/m value according to data
    v_e_div_m = float(sum(em_Values))/len(em_Values)
    v_slope_e_div_m = calcEm_wVSlope(v_slope, i_const)
    min = getBestEmValue(em_Values)
    best_v = [str(v_vs_R[min[0]][0])+" V", str(min[1]), str(min[3])]

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
                elif(i < 3):
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

def calcEm_wVSlope(b, i):
    return (125*(r_coil**2))/(32*(b**2)*(i**2)*(n**2)*(u**2))

def calcMagField(i):
    return (8*u*n*i)/(5*r_coil*sqrt(5))

def calc_e_div_m(v,b,r):
    return (2*v)/((b*r)**2)

def diffExpectedCalculated(calc):
    if(calc > e_div_m):
        return (calc - e_div_m)
    else:
        return (e_div_m - calc)

def getBestEmValue(em_Values):
    i = 0
    diff = diffExpectedCalculated(em_Values[i])
    min = [i, em_Values[i], diff, checkSig(diff)]
    for value in em_Values:
        diff = diffExpectedCalculated(value)
        if(diff < min[2]):
            min = [i, value, diff, checkSig(diff)]
        i+=1
    return min

def checkSig(diff):
    return (diff < e_div_m_unc)

def analyzeResult():
    print("With Varied Current, and Constant Voltage @ "+str(v_const)+" V:")
    printGraphCalc(i_slope, i_intercept, i_slope_e_div_m)
    printProgCalc(i_e_div_m, best_i)
    print("\n")
    print("With Varied Voltage, and Constant Current @ "+str(i_const)+" A:")
    printGraphCalc(v_slope, v_intercept, v_slope_e_div_m)
    printProgCalc(v_e_div_m, best_v)

def printGraphCalc(s, int, em):
    print(" Graphical Calculations:")
    print("     slope: "+str(s))
    print("     intercept: "+str(int))
    print("     slope-based e/m calculation: "+str(em))
    diff = diffExpectedCalculated(em)
    sig = checkSig(diff)
    print("     slope-based e/m value is significant: "+str(sig))

def printProgCalc(em, best_val):
    print(" Program Analysis:")
    print("     average e/m value from data: "+str(em))
    print("     best input for calculated e/m value: "+best_val[0])
    print("     e/m value for best input: "+best_val[1])
    print("     best e/m value is significant: "+best_val[2])

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
