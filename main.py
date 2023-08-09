import datetime

from chord_cutter import ChordCutter
from chord_cutter import segment_vertebrae


if __name__ == "__main__":
   a = datetime.datetime.now()

   cutter = ChordCutter("/path/to/dicom/ct",chord_name="chord_RS_name",seg_path="/path/for/TotalSegmentator/Outputs")

   cutter.add_vertebrae_group("cervical")    
   cutter.add_vertebrae_group("thorax")    

   segment_vertebrae(cutter.dcm_input_path,cutter.seg_path,cutter.vertebrae,fast=True)

   # To obtain only the slice index of the chord per vertebra
   # (no dicomRS is produced):
   # vdict = cutter.cut_chord(0.02,verbose=2)
   
   # To produce the dicomRS of the cutted chord
   cutter.vertebrae_RS_path = "/path/for/RTStruct/output"
   vdict = cutter.cut_chord(0.02,save_as_RS=True,verbose=2)

b = datetime.datetime.now() - a
print(b)