import argparse
import requests, sys
import pandas as pd
 

def main():
    server = "https://rest.ensembl.org"
    ext = "/cafe/genetree/member/id/ENSACAP00000022823"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        print('Reutrn not okay. Exit')
        sys.exit()
    
    decoded = r.json()
    keys = [k for k in decoded.keys()]
    print(keys)
    x = 1
    for t in decoded['tree']:
        print(t)
        print(decoded['tree'][t])

        if x == 4:
            break
        else:
            x += 1
        break
        



if __name__ == '__main__':
    main()