import plotly.graph_objs as go
import plotly.io as pio
from ast import literal_eval
from imageio import imread, mimsave
from os import remove
import numpy as np

def main():
	n = 1
	file_prefix = 'fig'
	color = 'Viridis'
	zvmax = float('-inf')
	zvmin = float('inf')
	images = []
	print('Converting output ...')
	with open("output_wave2d.txt", 'r') as f:
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
				
			fig = go.Figure()
			fig.add_heatmap(z=matrix, zmin=-32, zmax=20, colorscale=color)
			
			print(n)
			n += 1
			#pio.write_image(fig, f'{file_prefix}{int(counter + 1)}.png')
			images.append(pio.to_image(fig, format='png', width=700, height=500))
	
	print("Range: [", zvmin, ", ", zvmax, "]")
	print('Creating animation ...')
	convertGif(file_prefix, f'Wave2D_{color}', n, images)
	print('Cleaning up ...')
	cleanUp(file_prefix, n)
	print('Done.')
	
	
def convertGif(name_prefix, output_name, n, images):
	#images = []

	#for i in range(1,n):
	#	images.append(imread(f'{name_prefix}{i}.png'))
	mimsave(f'./{output_name}.gif', images)

def cleanUp(name_prefix, n):
	for i in range(1,n):
		remove(f"{name_prefix}{i}.png")
	
if __name__ == "__main__":
	main()


