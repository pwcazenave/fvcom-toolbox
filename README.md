#summary Basic Installation Instructions

=Dependencies =

You will need to download/install the netcdf toolbox which was originally written by Chuck Denham.  It still lives and functions with the Mathworks supplied mex-file interface (thanks to R. Signell for this info).  You can download the various netcdf interfaces using svn via the links below and adjust your matlab path accordingly.

{{{
  svn co http://mexcdf.svn.sourceforge.net/svnroot/mexcdf/mexnc/trunk
  svn co http://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/trunk
  svn co http://mexcdf.svn.sourceforge.net/svnroot/mexcdf/netcdf_toolbox/trunk
}}}

Make sure that "which mexnc"  returns an m-file, not a mex file.

=Optional Dependencies =

If you wish to project coordinates from geographic to Euclidean (see my_project.m) it is recommended you install the m_map matlab projection toolbox.
