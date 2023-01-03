def read_samplesCSV(spls):
	# reads in samples.csv file, format: Batch,Lane,Barcode,Sample,Alignments,ProviderName,Patient
	hdr_check = ['ExperimentFolder', 'Sample', 'AlignmentFolder', 'Path']
	switch = "on"
	file = open(spls, 'r')
	list_Experiment_Folder = []
	list_Sample = []
	list_AlignmentFolder = []
	list_Path = []
	for line in file:
		line = line.strip('\n').split(',')
		# Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
		if switch == "on":
			if (line == hdr_check):
				print("Passed CSV header check")
			else:
				Warning("CSV did NOT pass header check! Code continues, but first line ignored")
			switch = "off"
			continue
		# build lists
		list_Experiment_Folder.append(line[0])
		list_Sample.append(line[1])
		list_AlignmentFolder.append(line[2])
		list_Path.append(line[3])
	return [list_Experiment_Folder,list_Sample,list_AlignmentFolder,list_Path] # set(list_patient) provides only unique subject IDs