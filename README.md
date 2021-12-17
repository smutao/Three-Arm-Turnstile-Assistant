# Three-Arm Turnstile Assistant
PyMOL plugin to modify molecular structure via turnstile rotation




## Why we need this? 

We normally adjust molecular structure in terms of bond length, angle and dihedral. However, the present molecular structure editors (e.g. GaussView, Avogadro) cannot handle parameter like this 



<img src="https://user-images.githubusercontent.com/10540422/146474538-7cfdb858-bb4f-45f2-aeba-9d6721c50e5c.gif" alt="drawing" width="200"/>

which is called **three-arm turnstile rotation**. The fact is that such intramolecular motion is frequently observed in coordination compounds. We provide this tool *Three-Arm Turnstile Assistant* to help with the computational studies of this motion. 

## Installation 

1) clone this repo by clicking Code -> Download ZIP

2) open PyMOL 2.5+ 

3) click "Plugin" -> "plugin manager" -> "Install New Plugin" -> "Choose", then choose the "src/\_\_init\_\_.py" file

4) "Turnstile Assistant" will be installed to PyMOL and show up in the "Plugin" drop-in menu


## How to use? 

### Simple case (one atom per arm)


https://user-images.githubusercontent.com/10540422/146475323-773b76a4-fd32-41c0-b938-e4d5e0948c49.mp4





### Substituents on arm


https://user-images.githubusercontent.com/10540422/146475409-0fd4f3dc-9183-483f-9e75-f3972c742eee.mp4


**Note: For each arm, the atom directly bonded to the central anchor atom should be FIRST selected. In this video, the atoms directly bonded to the anchor atom are 31, 53, and 75.**



## Citation 

If you find this tool useful in your research, please cite our articles.






