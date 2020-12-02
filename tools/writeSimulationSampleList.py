
import argparse
import json

def main():
	
	args = getShellArguments()

	if (args.readsDir[-1] != "/"):
		args.readsDir = args.readsDir + "/"

	inputs = getSimulationInputs(args.inputJSON)

	writeSimulationSampleList(args, inputs)



def getShellArguments():
	'''
	Get shell enviroment arguments
	'''
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--inputJSON', 
		type=str, 
		help='ID of simulated cohort')
	parser.add_argument(
		'--readsDir',
		type=str,
		help='path to directory were the simulated reads files are saved')

	return parser.parse_args()


def getSimulationInputs(inputJSON):
	f = open(inputJSON,) 
	data = json.load(f) 	  
	f.close()
	return data

def writeSimulationSampleList(args, inputs):
	sampleListFileName = args.readsDir + inputs['cohort_id'] + ".samples.list"
	sampleList = open(sampleListFileName, "w")
	sampleList.write("# sample list for cohort: " + inputs['cohort_id'] + "\n")
	for i in range(inputs['n_samples']):
		sampleList.write(f"s_{i}\t{inputs['cohort_id']}.s_{i}.raw_1P.fq\t{inputs['cohort_id']}.s_{i}.raw_2P.fq\tILLUMINA\n")
	sampleList.close()


if __name__ == '__main__':
	main()