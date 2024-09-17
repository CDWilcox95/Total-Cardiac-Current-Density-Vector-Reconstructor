function fmdl=GetReconMesh(R,H, calc_new, cond_case, shape, Bh)

order=1;    L=32;
if calc_new
    cd eidors;
    run startup.m;  cd meshing/netgen;

    % Sphere specs:   (x0, y0, z0, r0)
%     sphere_specs=[Q0,0.1*R];  
%     heart_object_str="solid ball = sphere("+num2str(sphere_specs(1))+","+num2str(sphere_specs(2))+"," + num2str(sphere_specs(3))+";" +num2str(sphere_specs(4))+");";
%     extra={'ball',heart_object_str};

    
    
    %Cylinder Model Specs:  (cyl_shape, # electrodes, [# plane electrodes, ? , #
    %planes], electrode radius, extra interior specs)
    num_elec_plane=16;    num_rings=2;
    Z=[0.04, H-0.04];      

    ring_vert_pos=Z;   % Z planes where electrode rings are placed (cm)
%     elec_shape = [0.0254,0.0254,R*0.04];  
    elec_shape = [0.01, 0, R*0.04]; %Circular electrode configuration
    % elec_shape=[0,0,0.01];% point electrode configuration

    switch shape
        case "c"
            cyl_shape=[H, [R, [R*0.04]]];   % (height, [cylinder radius, max elem size])
            shift=2*pi/num_elec_plane;
            elec_pos=[(2*pi/num_elec_plane).*[1:num_elec_plane]'-shift, Z(1).*ones(num_elec_plane,1); (2*pi/num_elec_plane).*[1:num_elec_plane]'-shift, Z(2).*ones(num_elec_plane,1)];

            if strcmp(cond_case,"var")
                refine_cntr="("+num2str(Bh(1))+","+num2str(Bh(2))+","+num2str(Bh(3))+";"+num2str(Bh(4))+")";
                extra={'ball',"solid ball = sphere"+refine_cntr+" -maxh=0.003;"};
                fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape, extra);
            elseif strcmp(cond_case,"const")
                fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape);
            end

        case "e"
            ell_shape=[H, [0.15, 0.12, [R*0.04]]];
            shift=10*pi/num_elec_plane;
            elec_pos=[(2*pi/num_elec_plane).*[1:num_elec_plane]'-shift, Z(1).*ones(num_elec_plane,1); (2*pi/num_elec_plane).*[1:num_elec_plane]'-shift, Z(2).*ones(num_elec_plane,1)];

            [fmdl,mat_idx] = ng_mk_ellip_models(ell_shape, elec_pos, elec_shape);

    end



    % elec_pos = [16, Z];






    img= mk_image(fmdl,1); %img.elem_data(fmdl.mat_idx{2}) = 2;
    show_fem(img);

    cd ../../../;

    [~, vol_tetr]=CalculateMeshCoef(fmdl.elems, fmdl.nodes);


    if order==2
        fmdl.approx_type='tet10';
        [bound,elem,nodes] = fem_1st_to_higher_order( fmdl );   fmdl.boundary=bound;    fmdl.elems=elem;    fmdl.nodes=nodes;
        [fmdl.coef, ~]=CalculateMeshCoef(elem, nodes);
    end
    T=fmdl.elems;   Pt=fmdl.nodes;
    [C, ~]=CalculateMeshCoef(T,Pt);
    fmdl.nodes=Pt;

    
    fmdl.coef=C;    fmdl.volumes=vol_tetr;

    fmdl.area_elec=pi*elec_shape(1)^2;

    electrode_struct=fmdl.electrode;
    electrode_nodes=zeros(32,1);
    for l=1:32
        el=struct2cell(electrode_struct(l));   el_nodes=el{1,1};
        electrode_nodes(l)=el_nodes(1);
    end
    fmdl.electrode_idx=electrode_nodes;

    if order==2
        bndry_idx=fmdl.boundary;
        num_bndry_elems=length(fmdl.boundary);
        int_bndry=zeros(num_bndry_elems,1);
        for i=1:num_bndry_elems
            idx=bndry_idx(i,:);
            elem_pts=Pt(idx,:);

            for j=1:6
                if norm(elem_pts(j,:)-[0,0,H/2])<=0.2*H
                    int_bndry(i)=i;
                    break;
                end
            end

        end
        int_bndry=nonzeros(int_bndry);
        ext_bndry_idx=setdiff(1:num_bndry_elems, int_bndry);
        ext_bndry=bndry_idx(ext_bndry_idx,:);
        map=BndryElem2TetElem(ext_bndry, T);

        fmdl.face2elem=map;
    end
    fmdl.bndry2elecmap=GetElectrodeElems(fmdl);
    save("Cylinder_FEM_mesh_D"+shape+"_O"+num2str(order)+"_modeldata.mat", "fmdl");
else
    load("Cylinder_FEM_mesh_D"+shape+"_O"+num2str(order)+"_modeldata.mat"); 
    
    % cd eidors;
    % run startup.m;  cd meshing/netgen;
    % img= mk_image(fmdl,1); %img.elem_data(fmdl.mat_idx{2}) = 2;
    % show_fem(img);
end


end
