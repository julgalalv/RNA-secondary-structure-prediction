import os
import re
import sys
import argparse
import platform
import subprocess

'''
This program predicts the secondary structure of a RNA sequence
using CLINGO and maximizing two possible Energy functions. It also generates an image 
of the resulting structure using VARNA applet. Java is needed to run VARNA applet
More information about VARNA applet in http://varna.lri.fr/index.php?lang=en&page=home&css=varna

type in console 'python rna_prediction.py --help' to get the following info:

    usage: rna_prediction.py [-h] sequence energy

    Secondary structure prediction of a RNA sequence with image generation.

    positional arguments:
      sequence    RNA sequence [example: ACCGUA]
      energy      energy function used [possible values 1, 2]

    optional arguments:
      -h, --help  show this help message and exit

Based on 'Exploring Life through Logic Programming: Logic Programming in Bioinformatics -  RNA secondary 
structure prediction' available in https://computerscience.nmsu.edu/_files/documents/TR-CS-NMSU-2014-10-24.pdf
'''

### VARIABLES 
CLINGO_DIR = 'clingo'                                           # Clingo files directory
SYSTEM = platform.system().lower()                              # Executing system (windows, linux or darwin (MacOS))
CLINGO_SYS_DICT = {'windows':'clingo-5.4.0-win64', 
                   'linux' : 'clingo-5.4.0-linux-x86_64',
                   'darwin':'clingo-5.4.0-macos-x86_64'}

CLINGO_EXE_DIR = os.path.join(CLINGO_DIR,'bin',CLINGO_SYS_DICT[SYSTEM])  # Clingo executable directory
CLINGO_EXE = os.path.join(CLINGO_EXE_DIR,'clingo')                       # Clingo executable path

# Clingo lp file to be executed. The final value is modified using the console parameter 'energy'
ENERGY_FILE = 'rna_ss_prediction_E{}.lp'    

CLINGO_OUTPUT_DIR = os.path.join(CLINGO_DIR,'output')                   # Clingo output directory 
CLINGO_INPUT_DIR =  os.path.join(CLINGO_DIR,'input')                    # Clingo input directory
CLINGO_INPUT_FILE = os.path.join(CLINGO_INPUT_DIR,'rna_ss_input.lp')    # Clingo program input file (runtime generated)

# original jar from http://varna.lri.fr/index.php?lang=en&page=downloads&css=varna
VARNA_DIR = 'VARNA'                                             # VARNA file directory
VARNA_JAR = os.path.join(VARNA_DIR,'VARNAv3-93.jar')            # VARNA jar path
GENERATED_IMAGES_DIR = 'generated_images'                       # VARNA image generated directory

# Allowed bases.
BASES = 'ACGU'
CONSOLE_LINE_LENGTH_ = 80

### FUNCTIONS
def main():
    '''
    Main function. This program predicts the secondary structure of a RNA sequence
    using CLINGO and maximizing two possible Energy functions. It also generates an image 
    of the resulting structure using VARNA applet if Java is installed.
    '''
    global SEQUENCE     # Input sequence global variable
    global ENERGY_FILE  # ENERGY_FILE as global variable
    global ENERGY_FUNC  # Input energy function value global variable

    # Console argument parser
    parser = argparse.ArgumentParser(description='Secondary structure prediction of a RNA sequence with image generation.')
    parser.add_argument('sequence', type=str,
                         help='RNA sequence [example: ACCGUA]')
    parser.add_argument('energy', type=int,
                        help='energy function used [possible values 1, 2]')

    args = parser.parse_args() 

    SEQUENCE = args.sequence
    not_valid_seq = any(x not in BASES for x in SEQUENCE)
    if not_valid_seq: 
        print('Invalid sequence. Sequence can only have the bases A,C,G,U.')
        sys.exit()
    
    ENERGY_FUNC = args.energy
    if ENERGY_FUNC not in [1,2]: 
        print('Invalid energy function. Energy must be 1 or 2.')
        sys.exit()
    ENERGY_FILE = ENERGY_FILE.format(ENERGY_FUNC)

    print('='*CONSOLE_LINE_LENGTH_)
    print('SYSTEM PLATFORM USED: ',SYSTEM)
    print('CLINGO EXECUTABLE RELATIVE PATH: ',CLINGO_EXE)
    print('='*CONSOLE_LINE_LENGTH_,'\n')

    # Generates input file from input RNA sequence
    generate_input_file()

    # Executes Clingo program
    run_clingo()

    # Executes VARNA applet to generate image if java is installed
    if is_java_installed():
        generate_image()

def is_java_installed():
    '''
    Check if Java is installed to avoid executing VARNA jar if not.
    '''
    version_check = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT).decode("utf-8") 
    pattern = '\"(\d+\.\d+).*\"'
    version_search = re.search(pattern, version_check).groups()
    if len(version_search) != 0:
        version = version_search[0]
        print('USING JAVA VERSION: ',version)
        return True
    else:
        print('Java not found. Image will not be generated.')
        return False

def generate_input_file():
    '''
    Generates the file CLINGO_INPUT_FILE parsing a raw sequence into 
    clingo facts.

    Example: ACCG -> seq(1,a). seq(2.C). seq(3,c). seq(4,G)
    '''

    # File content instantation.
    content = ''

    print('PARSING SEQUENCE INTO CLINGO LP INPUT FILE...')
    for i,b in enumerate(SEQUENCE, start=1):
        new_line = '\n' if i % 6 == 0 else ''                       # line jump each 6 facts for better view
        content += 'seq({},{}). {}'.format(i,b.lower(),new_line)    # parsing
    # Write to file
    f = open(CLINGO_INPUT_FILE, 'w')
    f.write(content)
    f.close()
    print('Input File: ',CLINGO_INPUT_FILE,'\n')

def run_clingo():
    '''
    Run clingo program 'rna_ss_prediction_E{energy function value}.lp'. 
    The output is saved in the clingo output directory and showed in console.
    '''
    global CLINGO_OUTPUT_FILE
    # Output file name has the sequence at the end of the name 
    CLINGO_OUTPUT_FILE = ENERGY_FILE.split('.')[0]+'_{}.txt'.format(SEQUENCE)
    out_file = os.path.join(CLINGO_OUTPUT_DIR,CLINGO_OUTPUT_FILE)
    
    # Clingo program file path
    lp_file = os.path.join(CLINGO_DIR,ENERGY_FILE)
    
    # Run clingo program
    cmd = '{} {}'.format(CLINGO_EXE,lp_file,out_file)
    print('EXECUTING: ',cmd)
    print('RUNNING CLINGO...')
    clingo_output = subprocess.getoutput(cmd)
    
    # write output to file
    with open(out_file, "w") as clingo_output_file:
        clingo_output_file.write(clingo_output)
    clingo_output_file.close()
    
    # print output
    print('='*CONSOLE_LINE_LENGTH_)
    print(clingo_output,'\n')
    

def generate_image():
    '''
    Generate image from clingo output file using VARNA applet. To do so, we extract the pairing predicates
    from the clingo output file and construct the connection string needed by VARNA. This string has the length
    of the RNA sequence and a link is represented by parenthesis '(' for the first appearing base and ')'
    for the other one. Non paired bases are represented by dots '.'

    Example: 'ACCGA' with pairing(2,4) is represented by '.(.).'
    '''

    # Instance of the connection using dots
    connections = list('.'*len(SEQUENCE))
    # read Clingo output file
    with open(os.path.join(CLINGO_OUTPUT_DIR,CLINGO_OUTPUT_FILE),'r',encoding='utf-8') as clingo_out:
        content = clingo_out.read()
        # Get the number of models before optimum
        models = re.findall(r'(?<=Models       : )(.*)',content)
        # Extract the output between 'Answer: {models}' and 'OPTIMUM FOUND' (optimum model output)
        response = re.findall(r'(?<=Answer: {}\n)[\S\s]*(?=OPTIMUM FOUND)'.format(models[0]), content)
        # Extract pairing predicates from optimum model output
        pairings = re.findall(r'pairing\(\d+,\d+\)',response[0])
    clingo_out.close()

    # dictionary of the pairings (only once per pairing)
    pairing_dict = {}
    for p in pairings:
        indexes = re.findall(r'\d+',p) # Extract sequence indexes frfom pairing predicate
        n1 = int(indexes[0])
        n2 = int(indexes[1])
        if n1 < n2: pairing_dict[n1] = n2 # Add to dictionary 

    # Substitution in connection list
    for k,v in pairing_dict.items():
        connections[k-1] = '('
        connections[v-1] = ')'
    
    # connection to string
    connections = ''.join(connections)

    print('='*CONSOLE_LINE_LENGTH_)
    print('PAIRING DICTIONARY:      ', pairing_dict)
    print('VARNA SEQUENCE STRING:   ', SEQUENCE)
    print('VARNA CONNECTION STRING: ', connections)
    print('='*CONSOLE_LINE_LENGTH_, '\n')

    # print stats
    print_statistics(pairing_dict)

    # image output name and path
    image_name = SEQUENCE+'_'+str(ENERGY_FUNC)+'.png'
    image_file = os.path.join(GENERATED_IMAGES_DIR,image_name)
    image_title = "{}_E{}".format(SEQUENCE,ENERGY_FUNC)
    print('RUNNING VARNA. GENERATING IMAGE...')
    # VARNA command
    base_cmd = 'java -cp {} fr.orsay.lri.varna.applications.VARNAcmd'.format(VARNA_JAR)
    params = ' -sequenceDBN {} -structureDBN {} -title {} -o {}'.format(SEQUENCE,connections,image_title,image_file)
    # VARNA applet run
    os.system(base_cmd + params)

def print_statistics(pairing_dictionary):
    """
    Prints information about the sequence and the predicted
    secondary structure
    """
    sequence_dict = dict(enumerate(SEQUENCE, start=1))
    len_seq = len(SEQUENCE)
    pairings = len(pairing_dictionary)
    paired_bases = pairings * 2
    cg_pairings = 0
    au_pairings = 0
    gu_pairings = 0

    for k,v in pairing_dictionary.items():
        base_1 = sequence_dict[k]
        base_2 = sequence_dict[v]
        pairing = set([base_1,base_2])
        if pairing == set(['C','G']): cg_pairings += 1  
        if pairing == set(['A','U']): au_pairings += 1  
        if pairing == set(['G','U']): gu_pairings += 1  
    cg_prop = cg_pairings / pairings
    au_prop = au_pairings / pairings
    gu_prop = gu_pairings / pairings

    stats = {
        'SEQUENCE LENGTH:':[len_seq,'-'],
        'NUM. PAIRINGS:':[pairings,'-'],
        'NUM. PAIRED BASES:':[paired_bases,'-'],
        'NUM. CG PAIRINGS:':[cg_pairings,'-'],
        'NUM. AU PAIRINGS:':[au_pairings,'-'],
        'NUM. GU PAIRINGS:':[gu_pairings,'-'],
        'PROP. CG PAIRINGS:':[cg_prop, '0.53'],
        'PROP. AU PAIRINGS:':[au_prop, '0.35'],
        'PROP. GU PAIRINGS:':[gu_prop, '0.12'],

    }
    (print('STATS OF SEQUENCE: E{} - {}'.format(ENERGY_FUNC,SEQUENCE)))
    print('='*CONSOLE_LINE_LENGTH_)
    print ("{:<20} {:<20} {:<20}".format('', 'PREDICTED', 'EXPECTED'))
    for k, v in stats.items():
        pred, expected = v
        print ("{:<20} {:<20.2f} {:<20}".format(k, pred, expected))
    print('='*CONSOLE_LINE_LENGTH_,'\n')


if __name__ == '__main__':
    main()