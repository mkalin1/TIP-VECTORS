import numpy

import pandas


murphy_10={'LVIM':['LEU','VAL','ILE','MET'], 'ST':['SER','THR'] , 'FYW':['PHE','TYR','TRP'],'EDNQ':['GLU','ASP','ASN','GLN'], 'KR':['LYS','ARG']}


def main():
    #Important variables
    residue_stats_data={}
    tesselation_stats_data={}
    union_stats_data={}
    ratio_stats_data={}

    #1 P(Residue)
    for i in open('residue_stats_data.xls'):
        templine=i.strip().split()
        residue_stats_data[templine[0]]= float(templine[1])
    #Apply Murphy 10
    temp_residue_stats_data={}
    for i in residue_stats_data:
        flag=True
        for j in murphy_10:
            if(i in murphy_10[j]):
                try:
                    temp_residue_stats_data[j] += residue_stats_data[i]
                except:
                    temp_residue_stats_data[j] = residue_stats_data[i]
                flag=False 
        if(flag):
            temp_residue_stats_data[i] = residue_stats_data[i]
    residue_stats_data = temp_residue_stats_data
    #End Murphy 10
    residue_sum = numpy.sum( list(residue_stats_data.values()) )
    for i in residue_stats_data:
        residue_stats_data[i] = residue_stats_data[i] / residue_sum
    
    
    #2 P(Tesselations)
    for i in open('tesselation_stats_data.xls'):
        templine=i.strip().split()
        tesselation_stats_data[templine[0]]= float(templine[1])
    #Apply murphy 10
    temp_tesselation_stats_data = {}
    for i in tesselation_stats_data:
        
        new_key=[]
        for j in i.split('-'):
            flag=True
            for k in murphy_10:
                if(j in murphy_10[k]):
                    new_key.append(k)
                    flag=False
            if(flag):
                new_key.append(j)
        new_key_string='-'.join(sorted(new_key))
        try:
            temp_tesselation_stats_data[new_key_string] += tesselation_stats_data[i]
        except:
            temp_tesselation_stats_data[new_key_string] = tesselation_stats_data[i]
    tesselation_stats_data = temp_tesselation_stats_data
    #Murphy Ends
    tesselation_sum = numpy.sum( list(tesselation_stats_data.values()) )
    for i in tesselation_stats_data:
        tesselation_stats_data[i] = tesselation_stats_data[i] / tesselation_sum
    
    #3 P(Tesselations) / P(Residues)
    for i in sorted(tesselation_stats_data.keys()):
        union_stats_data[i] = numpy.prod([residue_stats_data[j] for j in i.split('-')])
        ratio_stats_data[i] = tesselation_stats_data[i] / union_stats_data[i]


    ranks = sorted(ratio_stats_data.items(), key=lambda kv: kv[1])

    #Plotting
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(95,95))
    
    #Ratio
    ax.plot( [i[1] for i in ranks] , [i[0] for i in ranks] , linewidth=2 ,linestyle='--', marker='o', color='b')
    
    
    ax.set(xlabel='P(Tesselation) / P(Res1).P(Res2).P(Res3).P(Res4)', ylabel='Tesselations')
    fig.tight_layout()
    plt.grid()
    fig.savefig('test.svg', format='svg', dpi=1200)
    #plt.savefig('test.eps', format='eps')
    plt.close()
    

    

    return True
if(__name__=='__main__'):
    main()