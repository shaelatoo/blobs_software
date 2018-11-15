# blobs_software

This is a partial collection of software I wrote in 2006-2008 for studying plasma blobs in COR1 images, dug up in 2018
to share.  I am including only a small portion of what is in my files that seems worthy of sharing.

I haven't tested or used most of it in a very long time.  Some dependencies may be missing.  Some of these files
may be obsolete versions whose more up-to-date versions I didn't find / have lost.  Caveat emptor.

The folder httm contains software used for creating height-time, moving difference images from raw COR1 fits files.
The folder radon constains software that uses the IDL radon transform to enhance height-time images.
The folder cluster contains higher-level analysis software.  The routine sector2.pro looks promising for taking 
  processed COR1 images (stored in a fits file?) and using the radon transform and my clustering software to locate
  and enhance upward-moving plasma blobs.
  
  
