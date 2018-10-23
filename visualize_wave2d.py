import plotly.graph_objs as go
import plotly.io as pio
from ast import literal_eval
from imageio import imread, mimsave
from os import remove
import numpy as np

def main():
	file_prefix = 'fig'
	color = 'Viridis'
	fname = "output_wave2d.txt"
	
	line_count = getLineCount(fname)
	print('Converting output ...', line_count)
	#matrices, zvmin, zvmax = processFile()
	print('Creating images ...')
	#createHeatmaps(matrices, color, zvmin, zvmax, file_prefix)
	print('Creating animation ...')
	#convertGif(file_prefix, f'Wave2D_{color}', len(matrices) + 1)
	print('Cleaning up ...')
	#cleanUp(file_prefix, len(matrices) + 1)
	print('Done.')


def getLineCount(fname):
	with open(fname, 'r') as f:
		for line_count, _ in enumerate(f):
			pass
		return line_count + 1
	
def processFile(fname):
	matrices = []
	zvmax = float('-inf')
	zvmin = float('inf')
	with open(fname, 'r') as f:
		lines = f.readlines()
		for counter, line in enumerate(lines[:-1]):
			line = line.replace('|\n', '')
			rows = line.split('|')
			matrix = []
			
			for row in rows:
				try:
					vals = list(literal_eval(row))
					matrix.append(vals)
					temp_min = np.min(vals)
					temp_max = np.max(vals)
					if(temp_min < zvmin):
						zvmin = temp_min
					if(temp_max > zvmax):
						zvmax = temp_max
					
				except Exception as e:
					print(row)
			
			matrices.append(matrix)
	return matrices, zvmin, zvmax

	
def createHeatmaps(matrices, color, zvmin, zvmax, file_prefix):
	print('In create Heatmaps, need to make ', len(matrices), ' heatmaps')
	for n, matrix in enumerate(matrices):
		print('Making figure ...')
		fig = go.Figure()
		print('Making new heatmap')
		fig.add_heatmap(z=matrix, zmin=zvmin, zmax=zvmax, colorscale=color)	
		print('Creating image')
		pio.write_image(fig, f'{file_prefix}{int(n + 1)}.png')
		print('Done with ', n)
	
	
def convertGif(name_prefix, output_name, n):
	images = []

	for i in range(1,n):
		images.append(imread(f'{name_prefix}{i}.png'))
	mimsave(f'./{output_name}.gif', images)

	
def cleanUp(name_prefix, n):
	for i in range(1,n):
		remove(f"{name_prefix}{i}.png")
	
if __name__ == "__main__":
	main()


