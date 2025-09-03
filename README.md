Title:Spatial patterns of trace fossil reveal escalating horizontal and vertical bioturbation in the terminal Ediacaran

Yarong Liu /yrliu@nigpas.ac.cn /yl2141@cam.ac.uk
State Key Laboratory of Palaeobiology and Stratigraphy, Nanjing Institute of Geology and Palaeontology, Chinese Academy of Sciences, Nanjing, 210008, China
University of Chinese Academy of Sciences, Beijing, 100049, China
Department of Zoology, and the University Museum of Zoology, University of Cambridge, Cambridge, CB2 3EJ, UK

Co-authors: Zhe Chen, Scott Evans, Shuhai Xiao, Xunlai Yuan, Emily G. Mitchell

Funding:
This research was supported by Science Fund for Creative Research Groups of National Natural Science Foundation of China (42130207)
This work has been supported by Natural Environment Research Council Independent Research Fellowship NE/S014756/1 to E.G.M. 
S.X. was supported by the National Science Foundation (EAR-2021207).
For the purpose of open access, the author has applied a Creative Commons Attribution 4.0 International (CC BY 4.0) license. licence to any Author Accepted Manuscript version arising. 
Citation:  
Licences: Published under a CC BY licence. 

#Data:
Contains data extracted from Ediacaran trace fossil Lamonte trevallis surface in the Shibantan biota, Yangtze Gorges area, South China. Length units are in centimeter. 

SData1 all tracefossils: data of traces, connected plugs and isolated plugs. 
	Containing positions(X, Y), axis length of plugs (Long axis, Short axis), orientation of plugs (Orientation), traces lengths (Length) and traces widths (Width). Extracted from Inkscape using a custom script written in Haskell (https://github.com/egmitchell/dex).

SData2 segment positions: data of traces' segments.
	Start (x0,y0) and end (x1,y1) points of all segments.
	Segments positions were extracted from Inkscape using a custom script written in Haskell (https://github.com/egmitchell/traces).  

Edge: positions of edge points

Data collection dates: Data extracted from photographs 2023. 

#Code
1.System requirements
	Operating system(s): Windows 11
	Programming language: R (version 4.3.2)
	Development environment: RStudio (version 2023.12 or later, optional but recommended)
	R packages:
		spatstat 3.1-1
		mclust 6.1.1
		methods 4.3.3
		plyr 1.8.9
		VGAM 1.1-11
		XML 3.99
		rlist 0.4.6.2
		ecespa 1.1-17
  	Hardware requirements:
  		Standard desktop or laptop computer
	
 2.Installation guide
 	Download data and code from https://github.com/yrLIU-ediacaran/Lamonte_sppa/tree/v1.0.1
  	Open R or RStudio (version 4.3.2 or later)
	Install required R packages in code, read csv file 
 	Installation time: less than 5 minutes on a standard desktop computer
  
3.demo
A minimal dataset (demo_dataset.csv) is provided for demonstration. This demo is a simplified subset of the full fossil dataset to allow quick testing of the code.
Running instruction:
	Open R or RStudio (version 4.3.2 or later)
	Load the required packages and demo script (demo.R)
Expected output
	Figure 1: Demo trace fossil point distribution. Circle sizes denotes the lengths of traces
	Figure 2: Density plot of points
	Figure 3: Pair correlation function (PCF) curve
	Figure 4: Mark correlation function (MCF) plot
 Expected run time
	~30–60 seconds on a standard desktop computer
