import sys
import phylotreelib as pt

#################################################################################################

def main():
    if len(sys.argv) != 7:
        print("usage: python parsemb.py trprobs parts vstat run1t run2t burninfrac")
        exit()
    trprobs_fname, parts_fname, vstat_fname, run1_fname, run2_fname, burninfrac = sys.argv[1:]
    burninfrac = float(burninfrac)
    id_leafdict = id_leafdict_from_trprobs(trprobs_fname)  # {id:leafname}
    id_branchdict = id_branchdict_from_vstat(vstat_fname)  # {id:[mean,var]}

    with open(parts_fname) as infile:
        partsdata = infile.readlines()
    for line in partsdata:
        if "*" in line:
            words = line.rstrip().split()
            parts_id = int(words[0])
            mean, var = id_branchdict[parts_id]
            astpos = asterisk_poslist(words[1])
            dotpos = dot_poslist(words[1])
            for bip1_id in astpos:
                leafname = id_leafdict[bip1_id]
                print(leafname, end=" ")
            print("|", end=" ")
            for bip2_id in dotpos:
                leafname = id_leafdict[bip2_id]
                print(leafname, end=" ")
            print("| {} {}".format(mean,var))

    tf1 = pt.Nexustreefile(run1_fname)
    tf2 = pt.Nexustreefile(run2_fname)
    trees1 = tf1.readtrees()
    trees2 = tf2.readtrees()
    discard1 = round(len(trees1)*burninfrac)
    discard2 = round(len(trees2)*burninfrac)
    trees1 = trees1[discard1:]
    trees2 = trees2[discard2:]
    with open("contest.postburnin.1.t", "w") as outfile:
        outfile.write(trees1.nexus())
    with open("contest.postburnin.2.t", "w") as outfile:
        outfile.write(trees2.nexus())

#################################################################################################

def id_leafdict_from_trprobs(trprobs_fname):
    """Parse mrbayes trprobs file to get translate block as a dict: {leafid:leafname}"""
    with open(trprobs_fname) as infile:
        data = infile.readlines()
    for i,line in enumerate(data):
        if "translate" in line:
            break
    leafdict = {}
    for j in range(i+1, len(data)):
        words = data[j].rstrip().split()
        if words[0] == "tree":
            break
        else:
            leafid = int(words[0])
            leafname = words[1][:-1]
            leafdict[leafid] = leafname
    return leafdict

#################################################################################################

def id_branchdict_from_vstat(vstat_fname):
    with open(vstat_fname) as infile:
        data = infile.readlines()
    id_branchdict = {}
    for line in data:
        if line.startswith("length"):
            line = line.replace("[", " ")
            line = line.replace("]", " ")
            words = line.rstrip().split()
            br_id = int(words[1])
            id_branchdict[br_id] = [float(words[2]), float(words[3])]  # [mean, var]
    return id_branchdict

#################################################################################################

def asterisk_poslist(asterisk_string):
    """Parse asterisk-dot string and return list of indices for asterisks"""
    indexlist = []
    for i,char in enumerate(asterisk_string):
        if char == "*":
            indexlist.append(i+1)
    return indexlist

#################################################################################################

def dot_poslist(asterisk_string):
    """Parse asterisk-dot string and return list of indices for asterisks"""
    indexlist = []
    for i,char in enumerate(asterisk_string):
        if char == ".":
            indexlist.append(i+1)
    return indexlist

#################################################################################################

if __name__ == "__main__":
    main()

