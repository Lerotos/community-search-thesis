#used to convert datasets obtained from https://snap.stanford.edu/data/index.html into a format useable by adj_list
#the hash function is used to turn the string of source and destination into an int
#the conversion doesn't seem to be perfect but it's good enough for now
#timestamp is shortened to just year for now
#was tested with the reddit dataset

with open("E:\\Uni\\bachelor\\1\\tempusL\\tempus\\t", "w") as writer:   #output file
    with open("E:\\Uni\\bachelor\\1\\tempusL\\tempus\\soc-redditHyperlinks-body.tsv", "r") as reader:   #input file
        for line in reader.readlines():
            templine = line.replace('\t', ' ')
            source = abs(hash(templine.split(' ')[0])//100000000000000)
            destination = abs(hash(templine.split(' ')[1])//100000000000000)
            time = templine.split(' ')[3].split("-")[0]
            #time = templine.split(' ')[3] + "T" + templine.split(' ')[4]
            writer.write("add " + str(source) + " " + str(destination) + " " + str(time) + "\n")
