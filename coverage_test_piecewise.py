
import pandas as pd
#import opt 
from sklearn.linear_model import LinearRegression
import numpy as np
import warnings
from getopt import getopt
import sys
import os

#option input
#option to change ranges of values 20-80

def Calculate_Slope(c,r):
    """Calculate Slope from a range.  Expects c to be the sample data array; and r to be a tuple for range which for calculate slope at."""
    X1 = np.array(range(r[0], r[1], 1)).reshape(-1, 1)
    Y = np.array(c.T.iloc[0, r[0]:r[1]])
    reg2 = LinearRegression().fit(X1, Y)
    m2 = reg2.coef_
    b2 = reg2.intercept_
    return (m2,b2)                
    
if __name__=='__main__':
    """Get summary info per position 20-80 for the 8 QC Samples.  Spec > 1.3 coverage may indicate bias"""
    #Defaults
    slope_spec_pass = 0.0055 #(*1.7=0.0046)
    slope_spec_fail = -0.005895
    section = (20, 80)
    
    #0-10
    lower_slope_spec_pass = (0.12,0.14) #5'
    input_coverage_file = 'compiled_matrix_5to3.txt'
    
    #Retrieve user parameters--> Here check how you can pass two values to r bug3pipeline???
    
    try:
        optlist, args = getopt(sys.argv[1:], "hi:r:")
    except:
        print("Error retrieving options")
        #print HELP_STRING
        sys.exit(1)

    for (opt, opt_arg) in optlist:
        if opt == "-h":
            #print ""
            #print HELP_STRING
            sys.exit(1)
        if opt == "-i":
            input_coverage_file = opt_arg
        if opt == "-r":
            #needs to be a tuple (start, end)
            section = tuple(opt_arg)


    qc = pd.read_table(input_coverage_file).drop('Normalized Gene Position', axis=1)
    sample_dict = dict()
    #print("Looking at files: %s" % (list(c.columns[1:])))
    transposed = qc.T

    #Calculate slope per sample(qc.columns)
    for c in qc.columns:
        sample = pd.DataFrame(qc[c])
        m,b = Calculate_Slope(sample, section)
        samplename = sample.columns[0]
        
        ##Create sections in coverage by looking at first 20 nt and last 20 nt
        lower_coverage = qc[c][0:20] 
        upper_coverage = qc[c][80:100]

        # Rationale is that, if the coverage slope in the fist 20 bases is greater than the slope of the uniform distribution we expect; a bias is present.
        # Also here we can pinpoint where the bias is present in the coverage length, and in what samples. 
        
        if (round(max(lower_coverage),2) > 1.4 or round(max(upper_coverage),2) > 1.4):
            ### max coverage or average of the coverage at 0-20 or 80-100 ?? I think average is better :/ 
            print(f"Sample: {samplename} failed with a slope of:{m} -5 or 3' bias'")
            passed= False

        case1 = (round(max(lower_coverage),2) >= 1.2 and round(max(lower_coverage),2) <= 1.4)
        case2 = (round(max(upper_coverage),2) >=1.2 and round(max(upper_coverage),2) <= 1.4 )
        
        elif case1 or case2 :
            m2, b2 = Calculate_Slope(sample, (0,20))
            m4, b4 = Calculate_Slope(sample, (80,100))
            #print(f"Slope of 0-10 section: {m2}")
            #print(f"Slope of 80-100 section: {m4}")
            
            #5'bias not enough to be detected
            if ( m2 <= lower_slope_spec_pass[1]):
                passed = True
                print(f"Sample: {samplename} passed with a slope of: {m}'")
            #3' bias not enough to be detected
            elif (abs(m4) <= lower_slope_spec_pass[1]):
                passed=True
                print(f"Sample: {samplename} passed with a slope of: {m}")
            else:
                passed = 'Need to confirm'
                print(f"Sample: {samplename} failed with a slope of: {m} -5 or 3' bias'")
            
        else: # less than 1.25
            if abs(m) > slope_spec_pass:
                print(f"Sample: {samplename} failed with a slope of: {m}")
                passed=False

            elif m < slope_spec_fail:
                print(f"Sample: {samplename} failed with a slope of: {m} -5 bias'")
                passed=False

            else:
                print(f"Sample: {samplename} passed with a slope of: {m}")
                passed = True
                m3, b3 = Calculate_Slope(sample, (0,10))
                m5, b5 = Calculate_Slope(sample, (80,100))
                #print(f"Slope of 0-10 section: {m3}")
                #print(f"Slope of 80-100 section: {m5}")
    
        sample_dict[c] = passed
        
    ## Add column to the all_metrics
    if 'all_metrics' in os.listdir('.'):
        all_metrics=pd.read_excel('all_metrics.xls')
        all_metrics.append(sample_dict,ignore_index=False)
    else:
        print(f'all_metrics not found in current directory {os.getcwd()}, not appending values')
    
sys.exit()                        
