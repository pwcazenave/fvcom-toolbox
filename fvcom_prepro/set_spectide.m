function set_spectide(Mobj,nComps,SpectralFile,MyTitle)

% Setup spectral tides on the open boundary and dump a spectral file
%
% function set_spectide(Mobj,nComps,SpectralFile,MyTitle)
%
% DESCRIPTION:
%    Setup spectral tides on the open boundary and dump a spectral file
%    This is a USER DEFINED driver program for the FVCOM spectral tide
%    It requires USER Modification to work
%
% INPUT
%    Mobj         = Matlab mesh object
%    nComps       = Number of tidal components
%    SpectralFile = Output file name
%    MyTitle     = Title in resulting NetCDF file.
%
% OUTPUT:
%
% EXAMPLE USAGE
%    set_spectide(Mobj,nComps,SpectralFile,MyTitle)
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2012-06-15 Added support for variables when calling set_spectide.
%	 2012-08-02 Can now write out equilibrium amplitudes and beta love numbers.
%
%==============================================================================
subname = 'set_spectide';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Set Component Periods
%------------------------------------------------------------------------------
% Components = {   'M2',    'N2',    'S2',   'O1',    'K1'};
% Period     = [44714.16, 45570.05, 43200, 92949.63, 86164.09];

%------------------------------------------------------------------------------
% Setup user defined phase and amplitude along the open boundaries
% need to set:
%   1.) Period - vector containing component period in seconds
%   2.) Amp    - array of size [Nobcs, Ncomponents] containing amplitude in m
%   3.) Phase  - array of size [Nobcs, Ncomponents] containing phase in degrees
%------------------------------------------------------------------------------
if nComps > numel(Mobj.Components)
    error('Too few components given in Mobj. Check and try again.')
end

if(Mobj.nObs==0)
	warning('cannot setup spectral open boundary, there is no open boundary in the mesh struct')
	return
end

cnt = 0;
Amp = nan(sum(Mobj.nObcNodes),nComps);
Phase = nan(sum(Mobj.nObcNodes),nComps);
ObcNodes = nan(1,sum(Mobj.nObcNodes));
for ob=1:Mobj.nObs
	nObcs = Mobj.nObcNodes(ob);
	for j=1:nObcs
		cnt = cnt + 1;
		ObcNodes(cnt) = Mobj.obc_nodes(ob,j);  % set the open boundary nodes

        Amp(cnt,:) = Mobj.amp_obc{ob}(:,j);
        Phase(cnt,:) = Mobj.phase_obc{ob}(:,j);
    end
end

%------------------------------------------------------------------------------
% Dump a spectral tide file in NetCDF
%------------------------------------------------------------------------------
%write_FVCOM_spectide(ObcNodes,Mobj.period_obc(1:nComps),Phase,Amp,Mobj.beta_love,Mobj.equilibrium_amp,SpectralFile,MyTitle)
write_FVCOM_spectide(ObcNodes,Mobj.Components,Mobj.period_obc(1:nComps),Phase,Amp,Mobj.beta_love,Mobj.equilibrium_amp,SpectralFile,MyTitle)
if(ftbverbose); fprintf(['end   : ' subname '\n']); end
