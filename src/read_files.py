NACA_6409_DIR = 'data\\NACA6409'
SG_6043_DIR = 'data\\SG6043'
FILES = ['50000', '100000', '200000', '500000', '1000000']
def main():
    output = open("data\\maximum_C_L_C_D.txt", "w")
    output.write("FILE     Re        alpha    C_L    C_D        CDp       CM     Top_Xtr    Bot_Xtr\n")
    ncrit = ''
    for _ in range(2):
        for f in FILES:
            NACA_path = NACA_6409_DIR + '\\xf-n6409-il-' + f + ncrit + '.txt'
            NACA_6409 = open(NACA_path, "r")
            best_line = ''
            ratio = -float('inf')
            index = 12
            for line in NACA_6409.readlines()[12: ]:
                l = line.split()
                '''alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr'''
                if (float(l[1]) / float(l[2])) > ratio:
                    ratio = float(l[1]) / float(l[2])
                    best_line = line
                index += 1
            output.write('NACA6409 ' + f + ncrit + best_line)
            NACA_6409.close()
        ncrit = '-n5'
    ncrit = ''
    for _ in range(2):
        for f in FILES:
            SG_path = SG_6043_DIR + '\\xf-sg6043-il-' + f + ncrit + '.txt'
            SG_6043 = open(SG_path, "r")
            best_line = ''
            ratio = -float('inf')
            index = 12
            for line in SG_6043.readlines()[12: ]:
                l = line.split()
                '''alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr'''
                if (float(l[1]) / float(l[2])) > ratio:
                    ratio = float(l[1]) / float(l[2])
                    best_line = line
                index += 1
            output.write('SG6043 ' + f + ncrit + best_line)
            SG_6043.close()
        ncrit = '-n5'
if __name__ == '__main__':
    main()