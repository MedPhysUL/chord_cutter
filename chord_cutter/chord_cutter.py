"""
   @file:               chord_cutter.py 
   @Author:             Louis Archambault
   
   @Creation Date:      08/2023 
   @Last modification:  08/2023
   
   @Description:        Contains class chord_cutter that manipulate dicom and
                        mask data to produce a segmented chord
"""

from pathlib import Path
from os.path import isdir, join, exists
from os import listdir
import tempfile

from rt_utils import RTStructBuilder
import numpy as np
import SimpleITK as sitk
import dicom2nifti
from totalsegmentator.python_api import totalsegmentator

class ChordCutter:
   """
   A class to mamage the preparation of sub-ROIs (one per vertebrae) of the
   spinal chord based on a DICOM CT dataset and DICOM-RT of the spinal chord
   """

   def __init__(
         self,
         dcm_input_path:str,
         ct_nifti_path:str = None,
         chord_name:str = None,
         seg_path:str = None 
   ):
      """
      Constructor of the ChordCutter class

      Parameters
      ----------
      dcm_input_path : str
         Path of the original dicom data (CT + RTStruct)
      chord_name : str
         Name of the chord structure in RTStruct
      ct_nifti_path : str
         Path of the nifti CT file generated from the dicom
      seg_path : str
         Folder where the segmented output will be placed
      """
      self._RS_path = None
      self._vertebrae_RS_path = None
      self.dcm_input_path = dcm_input_path
      if chord_name:
         self.chord_name = chord_name
      else:
         self._chord_name = "Moelle"
      self._ct_nifti_path = ct_nifti_path
      self._seg_path = seg_path
      self._vertebrae = []

   @property
   def dcm_input_path(self):
      return self._dcm_input_path

   @dcm_input_path.setter
   def dcm_input_path(self, dcm_path):
      if isdir(dcm_path):
         self._dcm_input_path = dcm_path
         for f in listdir(dcm_path):
            if f.startswith('RS'):
               self._RS_path = join(dcm_path,f)
      else:
         raise ValueError("invalid dcm input directory")
      
   @property
   def RS_path(self):
      return self._RS_path

   @property
   def chord_name(self):
      return self._chord_name

   @chord_name.setter
   def chord_name(self,name):
      if isinstance(name,str):
         self._chord_name = name
      else:
         raise TypeError("Chord name must be a string")
         
   @property
   def seg_path(self):
      return self._seg_path

   @seg_path.setter
   def seg_path(self, seg_path):
      if isdir(seg_path):
         self._seg_path = seg_path
      else:
         raise ValueError("invalid output directory")

   @property
   def vertebrae(self):
      return self._vertebrae

   def add_vertebrae(self,vert:str):
      if vert not in self._vertebrae:
         self._vertebrae.append(vert)
   
   def remove_all_vertebrae(self):
      self._vertebrae = []

   def add_vertebrae_group(self,vgroup:str):
      match vgroup:
         case "cervical":
            for v in ["vertebrae_C1","vertebrae_C2","vertebrae_C3","vertebrae_C4",
                      "vertebrae_C5","vertebrae_C6","vertebrae_C7"]:
               self.add_vertebrae(v)
         case "thorax":
            for v in ["vertebrae_T1","vertebrae_T2","vertebrae_T3","vertebrae_T4",
                      "vertebrae_T5","vertebrae_T6","vertebrae_T7","vertebrae_T8"
                      ,"vertebrae_T9","vertebrae_T10","vertebrae_T11","vertebrae_T12"]:
               self.add_vertebrae(v)
         case "lumbar":
            for v in ["vertebrae_L1","vertebrae_L2","vertebrae_L3","vertebrae_L4",
                      "vertebrae_L5"]:
               self.add_vertebrae(v)
         case _:
            raise ValueError("vertebrae groups: cervical, thorax, lumbar")
   
   @property
   def ct_nifti_path(self):
      return self._ct_nifti_path
   
   @ct_nifti_path.setter
   def ct_nifti_path(self,path):
      if path.endswith(".nii.gz"):
         self._ct_nifti_path = path
      else:
         raise ValueError("CT NIFTI path must end with .nii.gz")
   
   def get_ct_nifti(self) -> np.ndarray:
      if self._ct_nifti_path is None:
         raise ValueError("Must first set ct_nifti_path")
      
      create_niftii_from_dcm(self.ct_nifti_path,self.dcm_input_path)      
      return load_nifti(self._ct_nifti_path)

   @property
   def vertebrae_RS_path(self):
      return self._vertebrae_RS_path
   
   @vertebrae_RS_path.setter
   def vertebrae_RS_path(self,path:str):
      if isdir(path):
         self._vertebrae_RS_path = path
      else:
         raise ValueError("Must have a valid path for vertebrae RS")

   def cut_chord(self,T:float,save_as_RS=False, verbose=0) -> dict:
      """
         Function to divide the chord in segment based on the position of the
         vertebrae

         Parameters:
         -----------
         T : float
            A threshold representing the relative fraction of pixels in a slice
            relative to the total fraction of pixel in the vertebrae mask. Only
            slice with containing a fraction of pixels > T will be considered
         save_as_RS: bool
            If true, will produce dicomRS a file for each chord segment
         verbose : int
            Option to print some information as the function progresses
         
         Return:
         -------
         The output is a dict containing:
            {'<vertebrae name':(<idx of first CT slice>,<idx of last CT slice>)}
         If save_as_RS is True, dicom files will be saved in 
      """
      if save_as_RS:
         if self._vertebrae_RS_path is None:
            raise ValueError("Must define vertebrae_RS_path")
         chord = load_rs(self.dcm_input_path,self.RS_path,self.chord_name)
         all_seg_chord_masks = {}
         
      values = {}
      for v in self._vertebrae:
         full_path = join(self._seg_path,v+".nii.gz")
         if verbose > 1:
            print(full_path)
         vert = load_nifti(full_path,flip=True)

         slice = vert.sum((0,1)) # sum over, x, y
         total = vert.sum() # total pixels in mask
         if total <= 0:
            # This vertebrae wasn't segmented by TotalSegmentator
            pass
         else:
            arr_th = np.where((slice/total)>T)[0] # index of z coord where above T
            minv = arr_th.min()
            maxv = arr_th.max()
            values[v] = (minv,maxv)
            if verbose > 1:
               print(f'\t{v}, Min z: {minv}; Max z: {maxv}')
            if save_as_RS:
               chord_t = np.zeros(chord.shape)
               chord_t[:,:,minv:maxv] = chord[:,:,minv:maxv]
               if chord_t.sum() < 1:
                  if verbose > 0:
                     print(f'{v}: no overlap between chord and vertebrae')
               else:
                  all_seg_chord_masks[v] = chord_t.copy()
      if save_as_RS:
         vert_rtstruct = RTStructBuilder.create_new(dicom_series_path=self._dcm_input_path)
         for v in self.vertebrae:
            if v in all_seg_chord_masks.keys():
               msk = np.ma.getmask(np.ma.masked_greater(all_seg_chord_masks[v],0)) # convert to bool
               vert_rtstruct.add_roi(msk,name=v)
         vert_rtstruct.save(join(self._vertebrae_RS_path,"vertebrae.dcm"))
      return values

#---------------------------------------------------------------------------------
def segment_vertebrae(dcm_input_path,out_path,vertebrae_list,fast=False) -> bool:
   """
      A function that run TotalSegmentator based on information from ccutter
      return: True if TotalSegmentator was called, False otherwise
   """
   # get the segmentation that are already available
   seg_files = []
   for f in listdir(out_path):
      name = f.split('.')[0]
      seg_files.append(name)

   if all([v in seg_files for v in vertebrae_list]):
      print("All segmentations found, no need to run TotalSegmentator")
      return False
   else:
      print(f"Segmenting: {vertebrae_list}")
      totalsegmentator(input=dcm_input_path,output=out_path,output_type="nifti",fast=fast,roi_subset=vertebrae_list)
      return True

def create_niftii_from_dcm(nifti_path:str, input_dcm_path:str):
   """
      Uses dicom2nifti to produce a nifti image from a DICOM directory
      Supposed to convert a single image

      return False if no image was created
      return True if an image was created
   """
   if exists(nifti_path):
      print("nifti file exists, no conversion from DICOM needed")
      return False

   with tempfile.TemporaryDirectory() as tmpdirname:
      dicom2nifti.convert_directory(input_dcm_path, tmpdirname, compression=True, reorient=False)
      all_img_name = []
      for f in listdir(tmpdirname):
         if f.endswith(".nii.gz"):
            all_img_name.append(f)
      if len(all_img_name) > 1:
         print("More than one nifti image was produced, this shouldn't be. Only the first is used")
      elif len(all_img_name) < 1:
         raise FileExistsError("No nifti image created from dicom")
      else:
         ori_path = Path(join(tmpdirname,all_img_name[0]))
         destination = Path(nifti_path)
         ori_path.rename(destination)
   return 1

def load_nifti(path:str, flip=False) -> np.ndarray:
   """
      Load a nifti image (.nii.gz)

      Images produced by TotalSegmentator must be flipped This is not the case
      for those created by dicom2niftii without the reorient option
   """

   img = sitk.ReadImage(path)
   img_ar = sitk.GetArrayFromImage(img)

   # nifti images have 'z' as the first axe
   img_ar = np.swapaxes(img_ar,0,2)
   img_ar = np.swapaxes(img_ar,0,1)

   if flip:
      img_ar = np.flip(img_ar,0)
   
   return img_ar

def load_rs(dcm_path:str, rs_path:str, roi_name:str):
   rtstruct = RTStructBuilder.create_from(dcm_path,rs_path)
   mask = rtstruct.get_roi_mask_by_name(roi_name)
   return mask