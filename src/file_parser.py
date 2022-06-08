DIR = 'data\\SG6043\\xf-sg6043-il-100000.txt'
def main():
    out = open('data\\SG6043\\parsed.txt', 'w')
    with open(DIR, 'r') as file:
        data = file.readlines()
        l = len(data)
        for i in range(12, l):
            d = data[i]
            d = d.split('  ')
            out.write(d[1] + ' ' + d[2] + ' ' + d[3] + '\n')
    out.close()


if __name__ == '__main__':
    main()