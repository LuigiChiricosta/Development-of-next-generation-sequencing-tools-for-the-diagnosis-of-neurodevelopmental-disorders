def get_options():
        import optparse
        desc = 'Compact genome and exome GnomAD reduced datasets'
        parser = optparse.OptionParser("help: -f file1 -g file2", description=desc)
        parser.add_option("-f", "--file1", action="store", type='string', dest='f1', help='First file path')
        parser.add_option('-g', "--file2", action='store', type='string', dest='f2', help='Second file path')
        (options, args) = parser.parse_args()

	return args,options.f1,options.f2

def readFile(variantDictionary, filePath):
	f=0
	with open(filePath, "r") as fileName:
		for line in fileName:
			if f<0:
				return
			else:
				f += 1
				parameters = line.replace("\n", "").split("\t")
				ID = parameters[0].split(":")
				if len(ID) != 1:
					alternatives = ID[3]
					ACSplitted = parameters[1].split(",")
					HomSplitted = parameters[2].split(",")
					HemiSplitted = parameters[3].split(",")
					alternativesSplitted = alternatives.split(",")
					for a in range(len(alternativesSplitted)):
						alternative = alternativesSplitted[a]
						key = ID[0],ID[1],ID[2],alternative
						#print parameters, key
						AC = int(ACSplitted[a]) if a < len(ACSplitted) and ACSplitted[a] != "-" and ACSplitted[a] != "." else 0
						Hom = int(HomSplitted[a]) if a < len(HomSplitted) and HomSplitted[a] != "-" and HomSplitted[a] != "." else 0
						Hemi = int(HemiSplitted[a]) if a < len(HemiSplitted) and HemiSplitted[a] != "-" and HemiSplitted[a] != "." else 0
						if key not in variantDictionary:
							variantDictionary[key] = [AC, Hom, Hemi]
						else:
							variantDictionary[key] = [AC+variantDictionary[key][0], Hom+variantDictionary[key][1], Hemi+variantDictionary[key][2]]

if __name__ == '__main__':
	import sys

	args = sys.argv;
	args, file1, file2 = get_options()

	variantDictionary = {}
		
	if file1 != None and file2 != None:
		readFile(variantDictionary, file1)
		readFile(variantDictionary, file2)
	else:
		print "help: -f file1 -g file2	"

	print "ID\tAC\tHom\tHemi"
	for i in variantDictionary:
		print str(i[0])+":"+str(i[1])+":"+str(i[2])+":"+str(i[3])+"\t"+str(variantDictionary[i][0])+"\t"+str(variantDictionary[i][1])+"\t"+str(variantDictionary[i][2])+"\t"
