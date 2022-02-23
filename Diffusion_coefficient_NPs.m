HIHIHIHI HOSSEIN!
WHWJHosdfsdnfksjnflksfsdf,;lsd fldsk fljksnfjksnjfknskjflsdf
11111111111111111111111111111111111111111111111111111111111111
111111111111111111



0;
frame_number=2210;
box_max=1.3e2;
box_min=-2e0;
x_increment=.4;
length_box= box_max-box_min;
%     !!!!!!!!!!!!!!!     Format Should be like this length_box=1.3e2-(-2e0);
ratio_atoms_alighn=0.02;
atom_number_alighn=round(atom_number*ratio_atoms_alighn);

a = 2.755; %in A




textline_per_frame=9;
r_increment=.1;
tolerance_theta=2;
tolerance_x=1;
lattice_constant=3.52; % 2A is more than the true value in MD file
atom_radius=0.25*sqrt(2*(lattice_constant^2));  %just for FCC
probe_radius=2.5;
box_length=8.5826300000000003;

fid=fopen('AA_diffiusion_coefficient');
A=textscan(fid, '%s', 'delimiter', '\n');


%%%   Please NOTE that the MD file should be created sumetrically and the particle should be located at the origin of coordination system 




n=0;
for i=1:frame_number
    for j=1:atom_number
        
D=A{1}{n+textline_per_frame+j};
text_per_line=textscan(D, '%s', 'delimiter', ' ');

atom_coordinate(j,1,i)=str2double(text_per_line{1}{1});
atom_coordinate(j,2,i)=str2double(text_per_line{1}{2});
atom_coordinate(j,3,i)=str2double(text_per_line{1}{3});
atom_coordinate(j,4,i)=str2double(text_per_line{1}{4});
atom_coordinate(j,5,i)=str2double(text_per_line{1}{5});



    end
    n=n+atom_number+textline_per_frame;
end


mean_aggl=zeros(frame_number,3);

for i=1:frame_number
      for j=1:atom_number
          
mean_aggl(i,1)=mean_aggl(i,1)+atom_coordinate(j,3,i);
mean_aggl(i,2)=mean_aggl(i,2)+atom_coordinate(j,4,i);    
mean_aggl(i,3)=mean_aggl(i,3)+atom_coordinate(j,5,i);

      end
 
  
end
     mean_agglomerate=mean_aggl/atom_number;


for i=1:frame_number
      for j=1:atom_number

        find_id_1=find(atom_coordinate(:,1,i)==atom_coordinate(j,1,1))  ;      
atom_coordinate_trimmed(j,1,i)=atom_coordinate(find_id_1,1,i);
atom_coordinate_trimmed(j,2,i)=atom_coordinate(find_id_1,2,i);
atom_coordinate_trimmed(j,3,i)=atom_coordinate(find_id_1,3,i);
atom_coordinate_trimmed(j,4,i)=atom_coordinate(find_id_1,4,i);
atom_coordinate_trimmed(j,5,i)=atom_coordinate(find_id_1,5,i);


 end
      end
 
atom_coordinate=[];  
atom_coordinate=atom_coordinate_trimmed;
     
     
     

for i=2:frame_number
    
    
    
    for mno=1:atom_number
      dist_atoms(mno,i-1)=(((atom_coordinate(mno,3,i)-atom_coordinate(mno,3,1))^2)+((atom_coordinate(mno,4,i)-atom_coordinate(mno,4,1))^2)+((atom_coordinate(mno,5,i)-atom_coordinate(mno,5,1))^2));

    end
    sum_dist(1,i-1)=sum(dist_atoms(:,i-1));
    i
end
for i=2:frame_number
    
  

    diffusion_coeff_total(1,i-1)=sum_dist(1,i-1)/atom_number;
end






A=[];
atom_coordinate_trimmed=[];
atom_coordinate=[];

fid=fopen('AA_diffiusion_coefficient_total');
A=textscan(fid, '%s', 'delimiter', '\n');


n=0;
for i=1:frame_number
    for j=1:atom_number
        
D=A{1}{n+textline_per_frame+j};
text_per_line=textscan(D, '%s', 'delimiter', ' ');

atom_coordinate(j,1,i)=str2double(text_per_line{1}{1});
atom_coordinate(j,2,i)=str2double(text_per_line{1}{2});
atom_coordinate(j,3,i)=str2double(text_per_line{1}{3});
atom_coordinate(j,4,i)=str2double(text_per_line{1}{4});
atom_coordinate(j,5,i)=str2double(text_per_line{1}{5});



    end
    n=n+atom_number+textline_per_frame;
end


mean_aggl=zeros(frame_number,3);

for i=1:frame_number
      for j=1:atom_number
          
mean_aggl(i,1)=mean_aggl(i,1)+atom_coordinate(j,3,i);
mean_aggl(i,2)=mean_aggl(i,2)+atom_coordinate(j,4,i);    
mean_aggl(i,3)=mean_aggl(i,3)+atom_coordinate(j,5,i);

      end
 
  
end
     mean_agglomerate=mean_aggl/atom_number;


for i=1:frame_number
      for j=1:atom_number

        find_id_1=find(atom_coordinate(:,1,i)==atom_coordinate(j,1,1))  ;      
atom_coordinate_trimmed(j,1,i)=atom_coordinate(find_id_1,1,i);
atom_coordinate_trimmed(j,2,i)=atom_coordinate(find_id_1,2,i);
atom_coordinate_trimmed(j,3,i)=atom_coordinate(find_id_1,3,i);
atom_coordinate_trimmed(j,4,i)=atom_coordinate(find_id_1,4,i);
atom_coordinate_trimmed(j,5,i)=atom_coordinate(find_id_1,5,i);


 end
      end
 
atom_coordinate=[];  
atom_coordinate=atom_coordinate_trimmed;
     
     
     


for i=1:1
    
Atom=[];    
shp=[];    
polar_coord=[];
ijk=1;







 

% % % % % % 
% % % % % % centroid_x=0;
% % % % % % centroid_y=0;
% % % % % % centroid_z=0;
% % % % % % for mno=1:atom_number
% % % % % %   
% % % % % %     centroid_x=centroid_x+atom_coordinate(mno,3,i);
% % % % % %     centroid_y=centroid_y+atom_coordinate(mno,4,i);
% % % % % %     centroid_z=centroid_z+atom_coordinate(mno,5,i);
% % % % % % end
% % % % % %     centroid_x=centroid_x/atom_number;
% % % % % %     centroid_y=centroid_y/atom_number;
% % % % % %     centroid_z=centroid_z/atom_number;
% % % % % % 
% % % % % % for mno=1:atom_number
% % % % % %     
% % % % % %     atom_coordinate(mno,3,i)=atom_coordinate(mno,3,i)-centroid_x;
% % % % % %     atom_coordinate(mno,4,i)=atom_coordinate(mno,4,i)-centroid_y;
% % % % % %     atom_coordinate(mno,5,i)=atom_coordinate(mno,5,i)-centroid_z;
% % % % % % end

Atom=[atom_coordinate(:,3,i),atom_coordinate(:,4,i),atom_coordinate(:,5,i)];

AA(:,1,i)=Atom(:,1);
AA(:,2,i)=Atom(:,2);
AA(:,3,i)=Atom(:,3);

shp=[];
shp = alphaShape(Atom(:,1),Atom(:,2),Atom(:,3),3);


   for nn=1:atom_number

       
       
[polar_coord(1,nn),polar_coord(2,nn),polar_coord(3,nn)] = cart2sph(atom_coordinate(nn,3,i),atom_coordinate(nn,4,i),atom_coordinate(nn,5,i));

trial_movment_r=polar_coord(3,nn)+atom_radius;

[trial_movment_x,trial_movment_y,trial_movment_z]=sph2cart(polar_coord(1,nn),polar_coord(2,nn),trial_movment_r);




if inShape(shp,trial_movment_x,trial_movment_y,trial_movment_z)==0
    
   surface_atom_id(ijk,i)=atom_coordinate(nn,1,i); 
   surface_atom_row_id(ijk,i)=nn; 
    ijk=ijk+1;
    
end
   end







%disp('time steps=')
%i
end


for i=2:frame_number


    
    
    
    for mno=1:size(surface_atom_row_id,1)
 %  dist_atoms_surface(mno,i-1)=(((atom_coordinate(surface_atom_row_id(mno,1),3,i)-atom_coordinate(surface_atom_row_id(mno,1),3,1)^2)+((atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),4,i)-atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),4,1))^2)+(atom_coordinate(surface_atom_row_id(mno,1),5,i)-atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),5,1))^2)))));

     current_id=surface_atom_row_id(mno,1) ;
 dist_atoms_surface(mno,i-1)=(((atom_coordinate(current_id,3,i)-atom_coordinate(current_id,3,1))^2)+((atom_coordinate(current_id,4,i)-atom_coordinate(current_id,4,1))^2)+((atom_coordinate(current_id,5,i)-atom_coordinate(current_id,5,1))^2));
%        dist_atoms(mno,i-1)=(((atom_coordinate(mno,3,i)-atom_coordinate(mno,3,1))^2)+((atom_coordinate(mno,4,i)-atom_coordinate(mno,4,1))^2)+((atom_coordinate(mno,5,i)-atom_coordinate(mno,5,1))^2));

 
 %dist_atoms_surface(mno,i-1)=(((atom_coordinate(surface_atom_row_id(mno,1),3,i)-atom_coordinate(surface_atom_row_id(mno,1),3,1)^2)+((atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),4,i)-atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),4,1))^2)+(atom_coordinate(surface_atom_row_id(mno,1),5,i)-atom_coordinate(atom_coordinate(surface_atom_row_id(mno,1),5,1))^2)))));

    
      
    end
    sum_dist_surface(1,i-1)=sum(dist_atoms_surface(:,i-1));
    i

    



end




for i=2:frame_number
    
  
  diffusion_coeff_surface(1,i-1)=sum_dist_surface(1,i-1)/size(surface_atom_row_id,1);
end



% % % % % % % % 
% % % % % % % % 
% % % % % % % % fid = fopen('MD_1_trimmed5.txt', 'w'); 
% % % % % % % % 
% % % % % % % %   inner_atom_id=[];
% % % % % % % %    for imn=1:atom_number
% % % % % % % %     count_not_inner=0;
% % % % % % % %        
% % % % % % % %        for rst=1:size(surface_atom_id,1)
% % % % % % % %            
% % % % % % % %            if  atom_coordinate(imn,1,i)~=surface_atom_id(rst,i)
% % % % % % % %             count_not_inner=count_not_inner+1;   
% % % % % % % % 
% % % % % % % %        
% % % % % % % %            end
% % % % % % % %    
% % % % % % % %        
% % % % % % % %        end
% % % % % % % %          if count_not_inner==size(surface_atom_id(:,i),1)
% % % % % % % %               inner_atom_id=[inner_atom_id;atom_coordinate(imn,1,i)];
% % % % % % % %        
% % % % % % % %          end
% % % % % % % %    
% % % % % % % %    end
% % % % % % % %    
% % % % % % % %    
% % % % % % % %    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    for ij=1:size(surface_atom_id,1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         ABC=find(atom_coordinate==surface_atom_id(ij,i));   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        plot3(atom_coordinate(surface_atom_id(ij,i),3,i),atom_coordinate(surface_atom_id(ij,i),4,i),atom_coordinate(surface_atom_id(ij,i),5,i),'.','MarkerSize', 20,'color','b')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % fprintf(fid,"ITEM:"+" "+"TIMESTEP"+" "+"\n"+frame_number+"\n");
% % % % % % % % fprintf(fid,"ITEM:"+" "+"NUMBER OF ATOMS"+" "+"\n"+ atom_number+"\n");
% % % % % % % % fprintf(fid,"ITEM: BOX BOUNDS pp pp pp"+"\n");
% % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % fprintf(fid,-box_length+"e+01"+" "+box_length+"e+01"+"\n");
% % % % % % % % fprintf(fid,"ITEM:"+" "+"ATOMS"+" "+"id"+" "+"type"+" "+"x"+" "+"y"+" "+"z"+"\n");
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % for atom_id=1:atom_number
% % % % % % % % 
% % % % % % % %     number_not=0;
% % % % % % % %     for in=1: size(surface_atom_id,1)
% % % % % % % %         
% % % % % % % %     find_id=surface_atom_id(in,i)==atom_coordinate(atom_id,1,i)  ;    
% % % % % % % %     
% % % % % % % %      
% % % % % % % %      if find_id==1
% % % % % % % %              fprintf(fid,atom_id +" "+"1"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % %     
% % % % % % % %            
% % % % % % % %      end
% % % % % % % %      
% % % % % % % %  
% % % % % % % % 
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %      for in=1: size(inner_atom_id,1)
% % % % % % % %         
% % % % % % % %     find_id=inner_atom_id(in,1)==atom_coordinate(atom_id,1,i)  ;    
% % % % % % % %     
% % % % % % % %      
% % % % % % % %      if find_id==1
% % % % % % % %              fprintf(fid,atom_id +" "+"2"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % %     
% % % % % % % %            
% % % % % % % %      end
% % % % % % % %      
% % % % % % % %  
% % % % % % % % 
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % % % % % % %     for in=1: size(surface_atom_id,1)
% % % % % % % % % % % % %         
% % % % % % % % % % % % %     find_id=surface_atom_id(in,i)==atom_id  ;    
% % % % % % % % % % % % %     
% % % % % % % % % % % % %      
% % % % % % % % % % % % %      if find_id==0
% % % % % % % % % % % % %          number_not=number_not+1;
% % % % % % % % % % % % %              
% % % % % % % % % % % % %   
% % % % % % % % % % % % %       
% % % % % % % % % % % % %                 
% % % % % % % % % % % % %              
% % % % % % % % % % % % %      end
% % % % % % % % % % % % %      
% % % % % % % % % % % % %    if number_not== size(surface_atom_id,1) 
% % % % % % % % % % % % %      fprintf(fid,atom_id +" "+"2"+" "+atom_coordinate(atom_id,3,i)+" "+atom_coordinate(atom_id,4,i)+" "+atom_coordinate(atom_id,5,i)+"\n");
% % % % % % % % % % % % %    end
% % % % % % % % % % % % %      
% % % % % % % % % % % % %    
% % % % % % % % % % % % %  
% % % % % % % % % % % % % 
% % % % % % % % % % % % %     end
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % % 
% % % % % % % % end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 






% for i=1:size(diffusion_coeff_total,2)
%     
%   
% 
%     diffusion_coeff_total_time_ave(1,i)=sum(diffusion_coeff_total(1,1:i))/i;
% end

time_increment=100;

value_time_step=1e-3;

for i=1:frame_number-1

time_steps(1,i)=(i*time_increment*value_time_step);%%%Attention: time in picosecond!


end
plot(time_steps,diffusion_coeff_surface,'k', 'Linewidth', 2);

%xlabel('Time [ps]');
%xlim([.5e-9 1e-8]);
%ylabel('MSD [A^2]');
%ylim([1e-13 1e2])
%legend('Song et al(dp 3.52nm).','Song et al(dp 1.76nm).','This work(dp 3.52nm).','This work(dp 1.76nm).')

title('Ni/T=1000k')

