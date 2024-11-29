function [s,S2,CRV,COORD,NORMAL] = curv(Mask)
S2=[];CRV=[];COORD=[];NORMAL=[];
figure;
[ContN, hc]=contour(Mask,1);   % contour of the mask = to only get the border of the edge ?
% ContN = 2*n matrix whose two lines contain the x and y coordinates of the
% contour of the mask (i.e; plot(ContN(1,:), ContN(2,:)) sends the contour)
close(get( get( hc, 'Parent' ), 'Parent' ));
InGap=find(ContN(1,:)<1);   %find the x-coordinates of the contour  which are lower than 1 = gaps between the contour and the image border ?
NGap=size(InGap,2);
InGap(NGap+1)=size(ContN,2)+1;  %1st elem = number of gaps / 2nd one
for g=1:NGap
    Cont=ContN(:,InGap(g)+1:InGap(g+1)-2);  % = contour = plot(Cont(1,:), Cont(2,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Cont=(Cont(:,2:end-1));
    [Ds,Lin]=bwdist(~Mask);   %Ds = distance matrix of inside the spheroid = imshow(double(Ds/max(max(Ds))))
    [Ro,Co]=ind2sub(size(Mask),find(Ds==1));   %x,y coordinates of the edge  = scatter(Ro, Co,1)
    for n=1:size(Cont,2)   % for each pixel forming the contour
        [V,I]=min((Cont(1,n)-Co).^2+(Cont(2,n)-Ro).^2);   % we compute the distance from the other pixels and take the smallest one ?
        Cont(:,n)=[Co(I(1));Ro(I(1))];
    end
    % Front=bwperim(Mask);
    % Front([1 end],:)=0;Front(:,[1 end])=0;
    % [R,C]=find(Front);
    % Cont=[C';R'];
    [B,M,N]=unique(Cont','rows');
    Coord=Cont(:,sort(M));

     if size(Coord,2)>=9  %if the edge is delineated by more than 9 pixels
        Aux=[size(Coord,2)-3:size(Coord,2),1:size(Coord,2),1:4];
        for Cn=1:size(Coord,2)   % for each pixel ?
            D=abs(diff(Coord(1,Aux(Cn:Cn+8)))+i*diff(Coord(2,Aux(Cn:Cn+8))));%generates the length path
            s=[1,cumsum(D)+1];
            s2(Cn)=sum(D(2:3))/2;                                            %distance corresponding to the current point

            [Px,Sx]=polyfit(s,Coord(1,Aux(Cn:Cn+8)),2);
            [Py,Sy]=polyfit(s,Coord(2,Aux(Cn:Cn+8)),2);
            % Compute the derivative
            xder1 =2*Px(1)*s(5) + Px(2);
            yder1 =2*Py(1)*s(5) + Py(2);
            % construct the polynomial
            x=Px(1)*s(5)^2 + Px(2)*s(5) + Px(3);
            y=Py(1)*s(5)^2 + Py(2)*s(5) + Py(3);
            Crv(Cn)=(2*Py(1)*xder1-2*Px(1)*yder1)./(xder1^2+yder1^2).^1.5; % compute the radius of curvature
            Normal(Cn,:)=-[-yder1, xder1]/abs(-yder1+i*xder1); % Normalize the normal
             
        end
        S2=cat(2,S2,s2);
        CRV=cat(2,CRV,Crv);
        COORD=cat(2,COORD,Coord);
        NORMAL=cat(1,NORMAL,Normal);
    else
        disp('not enough points to compute the curvature');
        fprintf('\n');

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   clear Normal Crv s2
end

       
        
        
        
        
    
% % Coord=contour(Mask,1);
% % close;
% % Coord=(Coord(:,2:end-1));
% % if size(Coord,2)>=3
% %     Aux=[size(Coord,2),1:size(Coord,2),1];
% % 
% %     for Cn=1:size(Coord,2)
% %         a=Coord(:,Aux(Cn));
% %         b=Coord(:,Aux(Cn+2));
% %         c=Coord(:,Aux(Cn+1));
% %         v=b-a;
% %         if (dot(v,c-a)^2)/(dot(v,v)*dot(c-a,c-a)) == 1
% %             s(Cn)=sqrt(dot(v,v));
% %             s2(Cn)=s(Cn);
% %             Crv(Cn)=0;
% %         else
% %             lm=(dot(c,v)-dot(a,v))/dot(v,v);
% %             d=a+(lm*v);
% %             o=(2*d)-c;
% %             or(1:2,Cn)=o;
% %             r=sqrt(([a(1),c(1),b(1)]-o(1)).^2+([a(2),c(2),b(2)]-o(2)).^2);
% %             tet=atan2([a(2),c(2),b(2)]-o(2),[a(1),c(1),b(1)]-o(1));
% %             % % fit a second order parametric polynomial
% %             [Pr,Sr]=polyfit(tet,r,2);
% % 
% %             % Compute the derivative
% %             drdt =2*Pr(1)*tet(2) + Pr(2);
% %             ddrddt =2*Pr(1);
% %             r2=Pr(1)*tet(2)^2 + Pr(2)*tet(2) + Pr(3);
% %           
% %             s(Cn)=abs((Pr(1)/3)*tet(3)^3+(Pr(2)/2)*tet(3)^2+Pr(3)*(tet(3)-tet(1))-...
% %                 (Pr(1)/3)*tet(1)^3-(Pr(2)/2)*tet(1)^2);
% % s2(Cn)=sqrt(dot(c-a,c-a))+sqrt(dot(b-c,b-c));
% %             Crv(Cn)=(r2^2+(2*drdt^2)-(r2*ddrddt))/((r2^2+drdt^2)^(3/2))*(-1)^Mask(round(o(2)),round(o(1)));
% %          %        Crv(Cn)=abs(r2^2+(2*drdt^2)-(r2*ddrddt))/((r2^2+drdt^2)^(3/2));
% %         end
% %     end
% % end


% 
% for k=1:41%File.NFiles.PC
% 
%     
%      tempFileB = ['Boundary\automatic','\', 'manual_edgedist_crop_PC_', num2str(k), '.txt'];
%      EdgeDist=load(tempFileB);
%      Mask=logical(EdgeDist);
%       Cont=contour(Mask,1);
%          close;
%           Cont=(Cont(:,2:end-1));
%           [s,s2,Crv,or] = ultimate_curv(Mask);
%           ICv(k)=sum(Crv.*(s/2));
%           Per(k)=sum(s2/2);
%           SuC(k)=sum(Crv);
% end