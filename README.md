# Chord cutter

A simple tool that uses
[TotalSegmentator](https://github.com/wasserth/TotalSegmentator) to divide (or
*cut*) a spinal chord dicomRS structure into segments aligned on each vertebra.
The only input necessary are a CT dataset in DICOM and a RTStruct object of the
spinal chord. 

One possible workflow as used in `main.py` is as follow:
1. Run `TotalSegmentator` on the CT to segment each vertebra
2. Find the top and bottom CT slice for each vertebra
3. Find the interesection of the vertebral regions and the spinal chord
4. Save a RTStruct with the divided spinal chord

**Intended use**: for radiation therapy re-treatment, it is sometime necessary
to determine if there is a risk of overlap of dose hotspots within the spinal
chord. Because DVHs do not contain spatial information, dividing a spinal chord
ROI into sub-regions aligned with each vertebra, facilitate the evaluation of
possible dose overlaps.

Here is an example of the output for vertebrae C2 to T3 as visualized with [3D
Slicer](https://www.slicer.org/) (the amount of overlap between chord segments
can be controlled via the threshold parameter of `ChordCutter.cut_chord`):

![example of output](/imgs/chord_dcm.png)
