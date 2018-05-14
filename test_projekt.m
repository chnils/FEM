clear all
close all

[x,elem,nel,nnp,ndm,ndf,nen,matparam,drlt,robin,loadcurve,Q]=input_Rocket_Engine();

plotTriaMesh(elem,x);

% Alla Dofs, en för varje nod
allDofs = (1:1:(nnp*ndf))';

% dofs with Dirichlet boundary conditions
numDrltDofs = length(drlt.nodes);
drltDofs = zeros(numDrltDofs,1);
for i=1:numDrltDofs
    node = drlt.nodes(i);
    ldof = drlt.ldof;
    drltDofs(i) = (node-1)*ndf + ldof;
end

%FreeDofs
freeDofs = setdiff(allDofs, drltDofs);

% dofs with Robin boundary conditions
numRobinDofs = length([robin(1).nodes,robin(2).nodes,robin(3).nodes,robin(4).nodes,robin(5).nodes,]);
RobinDofs = zeros(numRobinDofs,1);
j = 1;
for k=1:5
    for i=1:length(robin(k).nodes)
        node = robin(k).nodes(i);
        ldof = 1;

        RobinDofs(j) = (node-1)*ndf + ldof;
        j = j + 1;
    end
end
%=========================================================================
% fe-analysis
%=========================================================================

% unknown node temperatures, begins at 20 Celcuis
a = 20*ones(ndf*nnp,1);

% %Temperature over time
% a_time=[];

% simulation time and time/load step
stopTime = 5000000.0;
dt = 100000;

% numerical tolerance to detect equilibrium
TOL = 1e-10;

% loop over all time/load steps
time = 0;
step = 0;
figure
while (time < stopTime)
    time = time + dt;
    % StopTime=1, => while-loopen körs 10 ggr.
    step = step + 1;   
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    
    % Robin values för robin noder. Uppdateras varje timestep. Oklart om
    % detta är rätt. Verkar inte användas till nåt.
    RobinValue = zeros(numRobinDofs,1);
    T_0 = matparam(5);
%     j = 1;
%     for k=1:5
%     for i=1:length(robin(k).nodes)
%         q_0 = 0;
%         %q_0 = (-1)*robin_1.alpha*(a(robin_1.nodes(i))-robin_1.Tinf);
%         %q_0 = (-1)*robin_1.alpha*(robin_1.Tinf);
%         RobinValue(j) = q_0*interp1(loadcurve(1).time, loadcurve(1).value, ...
%             time, 'linear', 'extrap');
%         j = j + 1;
%     end
%     end
    
    % Dirichlet boundary conditions
    drltValue = zeros(numDrltDofs,1);
    for i=1:numDrltDofs
        lcID = drlt.loadcurveID;
        scale = drlt.scale * loadcurve(lcID).scalefactor;
        drltValue(i) = scale*interp1(loadcurve(lcID).time, loadcurve(lcID).value, ...
            time, 'linear', 'extrap');
    end;
    
    %Intialise increasement in drlt
    a(drltDofs) = -150;
    drltValueIncr = drltValue - a(drltDofs);
    
    % always run through the while loop once to check equilibrium
    rsn = 1.0; % residuum norm
    iter = 0;
    
    while (rsn > TOL && iter < 2) %Denna loppen är skum, görs en gång
        %a 
        %Initialise global fint (internal flux)
        fint = zeros(ndf*nnp,1);
        
        %Initialise global fdyn (dynamic force vector)
        fdyn = zeros(ndf*nnp,1);
        
        %Initialise global fsur (surface vector)
        fsur = zeros(ndf*nnp,1);
        
        %Initialise global K (tangent stiffness matrix)
        K = zeros(ndf*nnp,ndf*nnp);
        
        %Initialise global C (Heat capacity matrix)
        C = zeros(ndf*nnp,ndf*nnp);
        
        %Initialise global rsd (residuum)
        rsd = zeros(ndf*nnp,1);
        
        
        
        %element loop
        for e=1:nel
            
            %dofs belonging to element e
            edof = elem(e).edof;
            
            % coordinates of the element nodes
            ex = x(elem(e).cn,1)*0.001;   % x-coord
            ey = x(elem(e).cn,2)*0.001;   % y-coord
            
            % area of the element 
            Ae=1/2*det([ones(3,1) ex ey]);
            
            % thickness of the plate 
            t = matparam(2);
            
            %Temperature gradient matrix,2D traingular element with 3 nodes
            B =  1/(2*Ae)*[ey(2)-ey(3)  ey(3)-ey(1)   ey(1)-ey(2); 
                           ex(3)-ex(2)  ex(1)-ex(3)   ex(2)-ex(1)];
            
            %Temperature values at the element nodes, to be updated in the
            %loop
            ae = a(edof);
            
            %Temperature values at element nodes from last iteration
            ae_prev = a(edof);
            
            % temperature gradient
            gradT = B*ae;
            
            % Fourier's law
            [q, D, elem(e).stateVar(1,:)] = fourierLaw2d_inl(gradT, matparam, elem(e).stateVar0(1,:));
            
       %----Skapa och fyll fsure som beror på robin villkoren-------------
          
            %Initialise local convection tangent stiffness matrix 
            Kce = zeros(3,3);
            
            %Initialise local fsur
            fsure = zeros(nen,1);
            
            BoolRobin = [0,0,0];
            for p=1:5
                if sum(ismember(elem(e).cn, robin(p).nodes))==2
                    BoolRobin = ismember(elem(e).cn,robin(p).nodes);
                    
                    if BoolRobin(1)==0
                        x1 = ex(2);
                        x2 = ex(3);
                        y1 = ey(2);
                        y2 = ey(3);
                    end
                     if BoolRobin(2)==0
                        x1 = ex(1);
                        x2 = ex(3);
                        y1 = ey(1);
                        y2 = ey(3);
                     end
                     if BoolRobin(3)==0
                        x1 = ex(1);
                        x2 = ex(2);
                        y1 = ey(1);
                        y2 = ey(2);
                     end
                     
                     %boundary length to be integrated over
                     L = sqrt((x2-x1)^2 + (y2-y1)^2);
                     
                      Kce_test = zeros(3,3);
                    if(BoolRobin == [1,1,0])
                        Kce_test = [2, 1, 0; 1, 2, 0; 0, 0, 0];
                    elseif(BoolRobin == [1,0,1])
                        Kce_test = [2,0,1;0,0,0;1,0,2];
                    elseif(BoolRobin == [0,1,1])
                        Kce_test = [0,0,0;0,2,1;0,1,2];
                    end
                    
                    Kce = L/6*Kce_test*robin(p).alpha*t;
                     
                    %Update local fsur vector
                    fsure = fsure + Kce*(ae-robin(p).Tinf);
                    line([x1 x2],[y1 y2]);
                    hold on;
                    
                end
            end
            
%             %index that will identify which boundary AE,EF,FI,IJ,JN that
%             %element is on
%             which_bound = 0;
%             double_boundaries = 0;
%             BoolRobin = [0,0,0];
%             k=1;
%             
%             %Här blir fel eftersom när vi kollar vilken boundary varje
%             %element tillhör så gör vi det i ordningen AE,EF,FI,IJ,JN
%             while(k<6)
%                 BoolRobin = ismember(elem(e).cn,robin(k).nodes);
%                 
%                 if sum(BoolRobin) == 2
%                     which_bound=k;
%                     k=6;
%                 end
%                 if sum(ismember(elem(e).cn,robin(k+1).nodes)) == 2
%                     double_boundaries = 1;
%                 end
%                 
%                 k=k+1;
%             end
%              
%             if which_bound~=0;       
%                     if BoolRobin(1)==0
%                         x1 = ex(2);
%                         x2 = ex(3);
%                         y1 = ey(2);
%                         y2 = ey(3);
%                     end
%                      if BoolRobin(2)==0
%                         x1 = ex(1);
%                         x2 = ex(3);
%                         y1 = ey(1);
%                         y2 = ey(3);
%                      end
%                      if BoolRobin(3)==0
%                         x1 = ex(1);
%                         x2 = ex(2);
%                         y1 = ey(1);
%                         y2 = ey(2);
%                      end
%                      
%                      %boundary length to be integrated over
%                      L = sqrt((x2-x1)^2 + (y2-y1)^2);
%                      
%                       Kce_test = zeros(3,3);
%                     if(BoolRobin == [1,1,0])
%                         Kce_test = [2, 1, 0; 1, 2, 0; 0, 0, 0];
%                     elseif(BoolRobin == [1,0,1])
%                         Kce_test = [2,0,1;0,0,0;1,0,2];
%                     elseif(BoolRobin == [0,1,1])
%                         Kce_test = [0,0,0;0,2,1;0,1,2];
%                     end
%                     
%                     Kce = L/6*Kce_test*robin(which_bound).alpha*t;
%                     
%                     %Update local fsur vector
%                     fsure = fsure + Kce*(ae-robin(which_bound).Tinf);
%                     line([x1 x2],[y1 y2]);
%                     hold on;
%             end
            
       %----fsure fylld----------------------------------------------------
       
       %----fyll i resterande lokala element vektorer/matriser--------------------
            % element C-matrix
            Ce = plantml(ex', ey', matparam(4)*matparam(3)); 
            %ska matparam(4) gångras med matparam(3)?
       
            %Local fint
            finte = -B'* q * Ae * t; %Kce*ae;
            
            % element dynamic force vector
            fdyne = (1/dt)*Ce*(ae-ae_prev)*t;
            
            % element stiffness matrix 
            Ke = (B')*D*B*t*Ae;
            
       %-------------------------------------------------------------------
       
       %----Assemble vectors och matrices into global variants------------
            
            % assemble finte into global internal force vector fint
            fint(edof) = fint(edof) + finte;
            
            % assemble fdyne into global dynamic force vector fdyn
            fdyn(edof) = fdyn(edof) + fdyne;
            
            % assemble fsure into global surface force vector fsur
            fsur(edof) = fsur(edof) + fsure;
                      
            % assemble Ke into the global stiffness matrix K
            K(edof,edof) = K(edof,edof) +Ke +Kce + Ce/dt;
           
        end
        %----end of element loop-------------------------------------------
    
        
        % calculate residuum and its norm
        rsd(freeDofs) = -fint(freeDofs) - fsur(freeDofs) - fdyn(freeDofs);           
        rsn = norm(rsd) + norm(drltValueIncr);
        fprintf(1, ' %2d. residuum norm= %e\n', iter, rsn);
        
        if (rsn >= TOL && iter < 1)
            disp('     solve and update');
            rsd(freeDofs) = rsd(freeDofs) - K(freeDofs,drltDofs)* drltValueIncr;
            % calculate increment of node displacements
            da = zeros(size(a));
            da(freeDofs) = K(freeDofs,freeDofs)\rsd(freeDofs);
            a(freeDofs) = a(freeDofs) + da(freeDofs);
            a(drltDofs) = a(drltDofs) + drltValueIncr;
            drltValueIncr = zeros(size(drltValueIncr));
        end
   
        iter = iter + 1;     
    end
    
    % save state variables for next time/load step
    for e = 1:nel
        elem(e).stateVar0 = elem(e).stateVar;
    end;
    
    postprocessing(step, time, ndm, elem, x, drltValue, drltDofs, fsur, RobinDofs, a);
       
end    


