function  postprocessingMovie(step, time, ndm, elem, X, drltValue, drltDofs, fsur, neumDofs, u)

% number of elements in the mesh
nel = length(elem);

%=========================================================================
% post-processing - visualising results
%=========================================================================

scale_factor = 1;

figure(1);
clf;
for e=1:nel
    nen = length(elem(e).cn);
    Xe = zeros(ndm, nen);
    for idm = 1:ndm
        Xe(idm,:) = X(elem(e).cn,idm);
    end
    
    % dofs belonging to element e
    edof = zeros(ndm, nen);
    for ien=1:nen
       for idm=1:ndm
           edof(idm,ien) = ndm*elem(e).cn(ien)+idm-ndm;
       end
    end
    ue = u(edof);
        
    xe = Xe + scale_factor * ue;
    
    % extract stress and strain tensor
    strain = elem(e).stateVar(1:4);
    stress = elem(e).stateVar(5:8);
    
    % plot deformed configuration  
    patch(xe(1,:), xe(2,:), stress(2));
    hold on

end
%
Xmax = max(X(:,1)); 
Xmin = min(X(:,1)) - 0.5;
Ymax = max(X(:,2)); 
Ymin = min(X(:,2));

title('sig_{yy}');
caxis('auto');
%h = colorbar;

axis equal;
hold off
drawnow;
pause(1);


