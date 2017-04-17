# Modified by Isaiah Solomon and Karol Josef Bustamante
# Original file name: stackdepthFasterNov15Version2.py 

# Originally Coded by Ananta Acharya
# Contact ananta@uga.edu if any questions
# This can be used as is or you can modify, The results are not guaranted to be correct, use with caution


import os, glob
from time import time
from scipy.misc import comb


def main():
    print "1. Make HAPLOTYPE and COUNT files from matches (For each SNP)"
    print "2. Define CRITERIA for FILTERING SNPs (From outputs of 1 above)"
    print "3. Change the FORMAT from A/C to binary for each SNP"
    print "4. Make HAPLOTYPE and COUNT files from matches (Haplotype format, like Stacks)"
    print "5. Filter the data"
    option = raw_input("Choose an OPTION: ")
    if option == "1":
        input_file              = raw_input("Where is the file-only DIRECTORY?\nAll matches file should be there: ")

        stackdepth(input_file)

    if option == "2":
        input_file_snp          = raw_input("Where is the SNP FILE processed from No. 1 option above?: ")
        input_file_count        = raw_input("Where is the COUNTS FILE processed from No. 1 option above?: ")
        minimum_homozygous      = raw_input("Minimum DEPTH to call homozygous?: ")
        minimum_heterozygous    = raw_input("Minimum DEPTH to call heterozygote (combined)?: ")
        heterozygous_ratio      = raw_input("Minimum RATIO to call heterozygous?: ")
        output_file             = raw_input("Output file NAME?: ")

        rehapmap(input_file_snp, input_file_count, output_file, minimum_homozygous, minimum_heterozygous, heterozygous_ratio)
        
    if option == "3":
        input_file              = raw_input("Where is the FILE processed from 1 or 2 above")
        output_file             = raw_input("Output file NAME?: ")

        changeformat(input_file, output_file)

    if option == "4":
        input_file              = raw_input("Where is the file-only DIRECTORY?\nAll matches file should be there: ")

        stackdepthHaplo(input_file)
    
    if option == "5":
        input_file              = raw_input("Where is the FILE processed from No. 2 option above?: ")
        output_file             = raw_input("Output file NAME?: ")
        filterMore(input_file, output_file)


def filterMore(input_file):
    ordered_dictionary     = collections.OrderedDict()
    output_file  = open(nameOfFile, "a")
    write_buffer = ""

    with open(input_file) as opened_file:

        for each_line in opened_file:
            line_list = each_line.strip().split(",")
            ordered_dictionary.update({line_list[0] : line_list})

    workspace = dict()

    for each_key in ordered_dictionary:
        key = ""

        for each_character in each_key:

            if each_character == '_':
                break

            key += character

        if not len(workspace):
            workspace[key] = each_key

        elif not (key in workspace) and len(workspace) == 1:
            key_, value_    = workspace.popitem()
            temp_list       = ordered_dictionary.get(value_)

            for index in range(len(temp_list)):
                write_buffer += temp_list[index] + " "

            write_buffer    += "\n"
            workspace[key]  = each_key

        elif key in workspace:
            workspace[key]  = each_key

        else:
            number_of_na    = 0
            na_list         = list()

            for eachKey in workspace:
                List = ordered_dictionary.get(workspace.get(eachKey))

                for each_word in List:

                    if (each_word == "NA"):
                        number_of_na += 1

                na_list.append(number_of_na)
                number_of_na = 0

            minimum = min(na_list)

            for index in len(na_list):

                if (na_list[index] == minimum):

                    for eachWord in workspace[index]:
                        write_buffer += eachWord + " "

                    write_buffer += "\n"
                    break

            workspace.clear()
            workspace[key] = each_key

    output_file.write(write_buffer)
    

def rehapmap(haplo_map, haplo_count, output_file, minimum_depth,minimum_depth_hetero, minimum_ratio):
    minimum_depth           = int(minimum_depth)-1
    minimum_depth_hetero    = int(minimum_depth_hetero)-1
    minimum_ratio           = float(minimum_ratio)
    all_pair_dictionary     = dict()
    haplo_file              = open(haplo_map,"r")
    rehap_file              = open(output_file, "a")

    for each_line in haplo_file:
        line_list                           = each_line.strip().split(",")
        all_pair_dictionary[line_list[0]]   = line_list

    with open(haplo_count) as haplo_count_file:
        line1           = next(haplo_count_file)
        write_buffer    = ""
        rehap_file.write(line1)

        for each_line in haplo_count_file:
            line_list       = each_line.strip().split(",")
            ID              = line_list[0]
            write_buffer    += ID

            for x in range(1, len(line_list) - 1):
                counts = line_list[x].split("|")

                if len(counts) == 1:

                    if counts[0] in ["NA"]:
                        to_write        = " NA"
                        write_buffer    += to_write
                        
                    else:
                        allele_depth = int(counts[0])
                        
                        if allele_depth > minimum_depth:
                            to_write        = all_pair_dictionary[ID][x]
                            write_buffer    += " " + to_write

                        else:
                            to_write        = " NA"
                            write_buffer    += to_write

                else:
                    allele1 = int(counts[0])
                    allele2 = int(counts[1])

                    if (allele1 + allele2) > minimum_depth_hetero:
                        snps = all_pair_dictionary[ID][x]
                        
                        if minimum_ratio <= (allele1 / float(allele1 + allele2)) <= (1 - minimum_ratio):
                            write_buffer += " " + snps

                        else:
                            snps_list       = snps.split("|")
                            minor_allele    = min(allele1, allele2)

                            if minor_allele == allele1:
                                minor_snp = snps[0]
                                major_snp = snps[1]

                            else:
                                minor_snp = snps[1]
                                major_snp = snps[0]
                            
                            p1 = binomial_test((allele1 + allele2), minor_allele)

                            if p1 >= 0.1:
                                write_buffer += " " + major_snp

                            elif p1<0.001:
                                write_buffer += " " + snps

                            else:
                                write_buffer += " NA"

                    else:
                        write_buffer += " NA"
                
            write_buffer += "\n"
        rehap_file.write(write_buffer)
    

def stackdepth(input_file):
    directory = input_file
    start = time()
    maxcat = 0

    haplo_dictionary    = dict()
    haplo_dictionary2   = dict()
    gencat_dictionary   = dict()
    snp_dictionary      = dict()

    os.chdir(directory)
    print "Looking in" +directory
    the_matching_files = glob.glob("*.matches.tsv")

    if len(the_matching_files) < 1:
        print "Match files not found"
        return

    gens = []

    for each_file in the_matching_files:
        print "parsing, " + each_file
        gen = each_file[:-12]
        gens.append(gen)

        the_file = open(each_file, "r")

        for each_line in the_file:
            words   = each_line.strip().split()
            catID   = words[2]
            maxcat  = max(maxcat, int(catID))
            haplo   = words[5]
            depth   = words[6]
            haplo_dictionary[gen, catID, haplo] = depth

            if not haplo.startswith("consensus"):

                for i in range(len(haplo)):
                    snp_dictionary[gen, catID, i, haplo[i]] = int(snp_dictionary.get((gen, catID, i, haplo[i]), 0)) + int(depth)

            gencat_dictionary[gen, catID] = 1

            if haplo not in haplo_dictionary2.setdefault(catID, []):
                haplo_dictionary2.setdefault(catID, []).append(haplo)

        the_file.close()

    writebuffer         = ""
    gens.sort()
    snp_output_file     = open(input_file + "/SNPAllmatches.snp.csv", "a")
    count_output_file   = open(input_file + "/SNPallmatches.count.csv", "a")
    snp_output_file.write("Catalog,")
    count_output_file.write("Catalog,")
    
    snp_output_file.write(",".join(gens))
    snp_output_file.write("\n")
    count_output_file.write(",".join(gens))
    count_output_file.write("\n")
    snp_buffer          = ""
    count_buffer        = ""
    count               = 0
    
    for cats in range(1, maxcat + 1):
        count   += 1
        haploo  = haplo_dictionary2.get(str(cats), 0) # print haploo
        
        if haploo<>0 and ("consensus" not in haploo):
            haploo = haploo[0]
        
            for i in range(len(haploo)):
                snp_linebuffer      = str(cats) + "_" + str(i) + ","
                count_linebuffer    = str(cats) + "_" + str(i) + ","

                for gen in gens:
                    snp_gen         = list()
                    count_gen       = list()
                    each_depth      = 0
                    eachgen_snpdict = dict()
                    flag            = gencat_dictionary.get((gen,str(cats)),0)

                    if flag == 1:
                        
                        for each_nucleotide in ["A","T","C","G"]:
                            depthh = snp_dictionary.get((gen,str(cats),i,each_nucleotide),0)

                            if depthh > 0:
                                eachgen_snpdict[each_nucleotide] = depthh
                                each_depth += depthh

                            if depth <= 0:
                                eachgen_snpdict[each_nucleotide] = 0

                        for key, value in eachgen_snpdict.items():
                            snp_gen.append(key)
                            count_gen.append(str(value))
                            
                        write_gensnp        = "|".join(snp_gen)
                        write_gencount      = "|".join(count_gen)
                        snp_linebuffer      += write_gensnp + ","
                        count_linebuffer    += write_gencount + ","

                    else:
                        snp_linebuffer      += "NA,"
                        count_linebuffer    += "NA,"
                        
                snp_buffer      += snp_linebuffer + "\n"
                count_buffer    += count_linebuffer + "\n"

            if count % 1000 == 0:
                current_time = time()
                new_time = str((current_time - start) / 60)
                
                print "%s rows in %s minutes" %(count, new_time)               
                    
                snp_output_file.write(snp_buffer)
                count_output_file.write(count_buffer)
                snp_buffer          = ""
                count_writebuffer   = ""

        
    current_time    = time()
    new_time        = str((current_time-start) / 60)
    
    print "%s rows in %s minutes" %(count, new_time)               
        
    snp_output_file.write(snp_buffer)
    count_output_file.write(count_writebuffer)
    snp_output_file.close()
    count_output_file.close()


def stackdepthHaplo(input_file):
    directory = input_file
    start = time()
    maxcat = 0

    haplo_dictionary   = dict()
    haplo_dictionary2  = dict()
    gencat_dictionary  = dict()

    os.chdir(directory)
    print "Looking in" + directory
    the_matching_files = glob.glob("*.matches.tsv")

    if len(the_matching_files) < 1:
        print "Match files not found"
        return

    gens = []

    for each_file in the_matching_files:
        print "parsing, " + each_file
        gen = each_file[:-12]
        gens.append(gen)
        the_file = open(each_file, "r")

        for each_line in the_file:
            words   = each_line.strip().split()
            catID   = words[2]
            maxcat  = max(maxcat, int(catID))
            haplo   = words[5]
            depth   = words[6]
            haplo_dictionary[gen, catID, haplo] = depth
            gencat_dictionary[gen, catID] = 1

            if haplo not in haplo_dictionary2.setdefault(catID, []):
                haplo_dictionary2.setdefault(catID, []).append(haplo)

        the_file.close()

    gens.sort()
    haplo_output_file   =open(directory + "/SNPAllmatches.haplo.csv","a")
    haplo_output_count  =open(directory +"/SNPallmatches.haplocount.csv","a")

    haplo_output_file.write("Catalog,")
    haplo_output_count.write("Catalog,")
    
    haplo_output_file.write(",".join(gens))
    haplo_output_file.write("\n")
    haplo_output_count.write(",".join(gens))
    haplo_output_count.write("\n")
    haplo_writebuffer   = ""
    count_writebuffer   = ""
    count = 0

    for cats in range(1,maxcat+1):
        count   += 1
        haploo  = haplo_dictionary2.get(str(cats), 0)

        if haploo<>0:

            if "consensus" in haploo:
                haplo_linebuffer    = str(cats) + ","
                count_linebuffer    = str(cats) + ","

                for gen in gens:
                    flag = gencat_dictionary.get((gen, str(cats)), 0)
                    
                    if flag == 1:  
                        depthh              = haplo_dictionary.get((gen, str(cats),"consensus"), 0)
                        write_genhaplo      = "consensus"
                        write_gencount      = str(depthh)
                        haplo_linebuffer    += write_genhaplo + ","
                        count_linebuffer    += write_gencount + ","

                    else:
                        haplo_linebuffer    += "NA,"
                        count_linebuffer    += "NA,"
                        
                haplo_writebuffer       += haplo_linebuffer + "\n"
                count_writebuffer       += count_linebuffer + "\n"
        
            else:
                haplo_linebuffer = str(cats) + ","
                count_linebuffer = str(cats) + ","

                for gen in gens:
                    haplo_gen           = list()
                    count_gen           = list()
                    eachgen_haplodict   = dict()
                    flag                = gencat_dictionary.get((gen, str(cats)), 0)

                    if flag == 1:

                        for haplotype in haploo:
                            depthh = haplo_dictionary.get((gen, str(cats), haplotype), 0)     
                            
                            if depthh > 0:
                                eachgen_haplodict[haplotype] = depthh

                        for key, val in eachgen_haplodict.items():
                            haplo_gen.append(key)
                            count_gen.append(str(val))
                                
                        write_genhaplo      = "|".join(haplo_gen)
                        write_gencount      = "|".join(count_gen)
                        haplo_linebuffer    += write_genhaplo + ","
                        count_linebuffer    += write_gencount + ","

                    else:
                        haplo_linebuffer    += "NA,"
                        count_linebuffer    += "NA,"
                        
                haplo_writebuffer   += haplo_linebuffer + "\n"
                count_writebuffer   += count_linebuffer + "\n"

            if count % 1000 == 0:
                current    = time()
                new_time   = str((current - start) / 60)
                
                print "%s rows in %s minutes" % (count, new_time)               
                    
                haplo_output_file.write(haplo_writebuffer)
                haplo_output_count.write(count_writebuffer)
                haplo_writebuffer   = ""
                count_writebuffer   = ""

    current    = time()
    new_time   = str((current - start) / 60)
    
    print "%s rows in %s minutes" % (count, new_time)               
        
    haplo_output_file.write(haplo_writebuffer)
    haplo_output_count.write(count_writebuffer)
    haplo_output_file.close()
    haplo_output_count.close()


def changeformat(input_file, output_file):
    """Changes format from SNP A/C...to binary"""
    inputed_file    = open(input_file, "r")
    header          = next(inputed_file)
    write_buffer     = header
 
    
    for each_line in inputed_file:
        line_list   = each_line.strip().split(",")
        catalog     = line_list[0]
        snps        = line_list[1::]

        for snp_gen in snps:

            if "|" in snp_gen:
                nuclei = snp_gen.split("|")
                break

        for nucleus in nuclei:
            each_line = catalog + "_" + nucleus
            
            for snp_gen in snps:

                if "NA" in snp_gen:
                    each_line += ",NA"

                elif nucleus in snp_gen:
                    each_line += ",1"

                else:
                    each_line += ",0"

            write_buffer += each_line + "\n"

    print "writing in %s.f01.csv" % inf
    opened_file = open(output_file,"a")
    opened_file.write(write_buffer)
    inputed_file.close()
    opened_file.close()


def binomial_test(n, k):
    p1 = comb(n, k) * 0.0025 ** k * 0.9975 ** (n - k)
    return p1

main()
