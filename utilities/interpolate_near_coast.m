function [data]=interpolate_near_coast(distbc,tempzet,doExtrap)
% we should be using the law of the wall here s
% to interpolate the velocitie
        [~,imaxd]=max(distbc);
        tempzet1=tempzet(1:imaxd);
        distbc1=distbc(1:imaxd);
        % find nans and non nans
        ibadz=find(isnan(tempzet1));
        igoodz=find(~isnan(tempzet1));
        if doExtrap
        % Uses csaps as a smoothing spline interpolant
        pp=csaps(distbc1(igoodz),tempzet1(igoodz));
        % fnxtr allows usage of interpolan spline for extrapolation
        pp = fnxtr(pp);
        tempzet1(ibadz)=fnval(pp,distbc1(ibadz));
        else
            % use good point closest to coast
            [~,use_idx]=min(distbc1(igoodz));
           tempzet1(ibadz)= tempzet1(igoodz(use_idx));
        end
        % evaluate spline on nan points near the coast.
        % Repeat for the second stretch of boundary
        tempzet2=tempzet(imaxd+1:end);
        distbc2=distbc(imaxd+1:end);
        ibadz=find(isnan(tempzet2));
        igoodz=find(~isnan(tempzet2));
                if doExtrap

        pp=csaps(distbc2(igoodz),tempzet2(igoodz));
        pp = fnxtr(pp);
        tempzet2(ibadz)=fnval(pp,distbc2(ibadz));
        else
            % use good point closest to coast
            [~,use_idx]=min(distbc2(igoodz));
           tempzet2(ibadz)= tempzet2(igoodz(use_idx));
        end
        data=[tempzet1(:);tempzet2(:)];
return
