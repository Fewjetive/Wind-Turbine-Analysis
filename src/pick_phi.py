import pandas as pd
import numpy as np
from math import pi

def main():
    directory = "data\\results\\" 
    filename = "100x1000_14"
    excel_name = directory + filename + '.xlsx'
    phi = pd.read_excel(excel_name, header=None, sheet_name=0)
    a = pd.read_excel(excel_name, header=None, sheet_name=1)
    b = pd.read_excel(excel_name, header=None, sheet_name=2)
    phi = phi.to_numpy(dtype=float)
    a = a.to_numpy(dtype=float)
    b = b.to_numpy(dtype=float)
    (row, col) = phi.shape
    print(phi.shape)
    np_arr = np.zeros((row, 4), dtype=float)
    i = 0
    for r in range(row):
        for c in range(col):
            if c < col - 1 and not np.isnan(phi[r][c+1]):
                continue
            np_arr[i][0] = phi[r][c]
            np_arr[i][1] = a[r][c]
            np_arr[i][2] = b[r][c]
            np_arr[i][3] = (2 * (pi * ((r + 1) / row)) * ((c + 1) / 2 / col)) / 3
            i += 1
            break
    df = pd.DataFrame(np_arr, columns=['phi', 'a', 'b', 'c'])
    with pd.ExcelWriter(directory + 'out_' + filename + '.xlsx', engine='xlsxwriter') as writer:
        df.to_excel(writer)
if __name__ == '__main__':
    main()