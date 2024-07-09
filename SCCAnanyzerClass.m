classdef SCCAnanyzerClass
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Identify Neighbors and Neighbor Structure

        %identify the neighbor structure
        % omega = depth of the in-neighbor structure
        % disc = whether the graph is strongly connected
        % xi(:,:,n+1) = connectivity from in-neighbors
        function [omega,disc,xi]=NeighborStructure(G,n)
            disc=true;
            xi(:,:,1)=eye(n);
            xit=eye(n);
            xip=eye(n);

            for t=2:1:n+1
                for i=1:1:n
                    N=nearest(G,i, 1, "Direction","incoming");
                    for j=1:1:n
                        for k=1:1:n
                            if (ismember(k,N)>0)
                                xit(i,j)=max(xit(i,j),xip(k,j));
                            end
                        end
                    end
                end
                xi(:,:,t)=xit;
                xip=xit;
            end

            omega(:,:,1)=zeros(n,n);

            for t=2:1:n+1
                omega(:,:,t)=omega(:,:,t-1)+(t-1)*(xi(:,:,t)-xi(:,:,t-1));
            end

            for i=1:1:n
                for j=1:1:n
                    if (i~=j) && (omega(i,j,n+1)==0)
                        disc=false;
                    end

                end
            end

        end

        %% Convert into Latex Format for easy Printing
        %-------------------------------
        %Write a matrix into LaTeX file
        function s=MatrixOutput(W)
            
            % Convert
            % Get matrix dimensions
            m = size(W, 1);
            n = size(W, 2);
            % Create first line
            s = sprintf('  \\begin{bmatrix}\n  ');
            % Add matrix content
            for k = 1:m
                for l = 1:n
                    s = sprintf('%s %6.0f', s, W(k, l));
                    if l < n
                        s = sprintf('%s &', s);
                    end
                end
                if k < m
                    s = sprintf('%s \\cr', s);
                end
                s = sprintf('%s\n  ', s);
            end
            % Add last line
            s = sprintf('%s\\end{bmatrix}\n', s);
            % Print the result if needed
            % disp(s);
        end

        %% Find information number


        %-----------------------------------
        %identify the information number
        %Information number is the sum of known neighbors +1(itself)
        % ci(:,:,n+1) = rows of \eta_i(n-1) in the slide
        function [ci]= InformationNumber(G,n)
            [omega,disc,xi] = SCCAnanyzerClass.NeighborStructure(G, n);
            chi=sum(xi(:,:,n+1),2);
            ci = eye(n) .* chi;
            cit=ci;
            cip=ci;
            for t=2:1:n+1
                for i=1:1:n
                    N=nearest(G,i, 1, "Direction","incoming");
                    for j=1:1:n
                        for k=1:1:n
                            if (ismember(k,N)>0)
                                cit(i,j)=max(cit(i,j),cip(k,j));
                            end
                        end
                    end
                end
                ci(:,:,t)=cit;
                cip=cit;
            end
        end

        % ----------------------
        %% Find The SCCs by identifying information number
        %Find the SCCs
        % SC = list of SCCs
        % SP = list of in-neighbors for each node and its SCC
        function [SC, SP]=findingSCC(G,n)
            [ci]= SCCAnanyzerClass.InformationNumber(G,n);
            ci=ci(:,:,n);
            SC = cell(1, n);
            SP = cell(1, n);
            for i = 1:n
                sci = [];
                spi=[];
                for j = 1:n
                    if ci(i, j) == ci(i, i)
                        sci = [sci, j];
                    end
                    if 0 < ci(i, j) && ci(i, j) < ci(i, i)
                        spi = [spi, j];
                    end
                end
                SC{i} = sci;
                SP{i} = spi;
            end
        end

        %% Find if the SCC is a sink or source
        %sink = information comes in (but doesn't go out)
        %Source = information goes out (but doesn't come in)
        

        %-------
        %Find if the SCC is a source or sink
        function [oi]= findingsource(G,n)
            oi=zeros(n);
            [SC, SP]=SCCAnanyzerClass.findingSCC(G,n);
            for i=1:1:n
                No=nearest(G,i, 1, "Direction","outgoing");
                for j=1:1:n
                    if i==j
                        for k=1:1:n
                            if (ismember(k,No)>0) && ~(ismember(k,SC{i})>0)
                                oi(i,i,1)=1;
                            end
                        end
                    end
                end
            end

            oit=oi(:, :, 1);
            oip=oi(:, :, 1);
            for t=2:1:n+1
                for i=1:1:n
                    N=nearest(G,i, 1, "Direction","incoming");
                    for j=1:1:n
                        for k=1:1:n
                            if (ismember(k,N)>0)
                                oit(i,j)=max(oit(i,j),oip(k,j));
                            end
                        end
                    end
                end
                oi(:,:,t)=oit;
                oip=oit;
            end
        end


        %---------------------
        %% Index accoding to SCC nature

        %Index and group the representative nodes
        % Vsour = indices of source SCCs
        % Vsink = indices of sink SCCs
        % Visol = indices of isolated SCCs
        % Vmid = indices of SCCs that have both in-neighbor SCCs and
        %           outneighbor SCCs
        function [Vsour, Vsink, Visol,Vmid] = indexSCC(G,n)
            [SC, SP]=SCCAnanyzerClass.findingSCC(G,n);
            [oi]= SCCAnanyzerClass.findingsource(G,n);
            oi=oi(:,:,n);
            vso=[];
            vsi=[];
            vis=[];
            vsu=[];
            flg=false;
            flg1=false;
            for i=1:1:n
                mxSC=max(SC{i});
                for j=1:1:n
                    if (ismember(j, SC{i})) && (oi(i,j)==1)
                        flg=true;
                        break;
                    else
                        flg=false;
                    end
                end
                if flg==true && all(SP{i}==0)
                    vso=[vso mxSC];
                end
                if flg==false && ~all(SP{i}==0)
                    vsi=[vsi mxSC];
                end
                if flg==false && all(SP{i}==0)
                    vis=[vis mxSC];
                end
                if flg==true && ~all(SP{i}==0)
                    vsu=[vsu mxSC];
                end

            end
            Vsour=unique(vso);
            Vsink=unique(vsi);
            Visol=unique(vis);
            Vmid=unique(vsu);
        end


        %------------------------------

        
        %% Find SCC Structure

        % Analyze SCC structure
        % theta = in-neighbor structure for SCCs
        % eta = depth of in-neighbor structure for SCCs
        function [theta, eta] = SCCStructure(G, n)
            theta(:,:,1)= zeros(n, n);
            eta(:,:,1) = zeros(n, n);

            [SC, ~] = SCCAnanyzerClass.findingSCC(G, n);

            % Initialize theta matrix
            for i = 1:n
                for j = 1:n
                    if ismember(j, SC{i})
                        theta(i, j, 1) = 1;
                    end
                end
            end

            theit = theta(:,:,1);
            theip = theta(:,:,1);
            for t = 2: n+1
                for i = 1:n
                    N = nearest(G, i, 1, 'Direction', 'incoming');
                    for j = 1:n
                        for k = 1:n
                            if (ismember(k, N)>0)
                                theit(i, j) = max(theit(i, j), theip(k, j));
                            end
                        end
                    end
                end
                theta(:, :, t) = theit;
                theip = theit;
            end

            % Initialize eta matrix
            etaz = theta(:, :, 1);
            etax=theta(:, :, 1);
            eta(:, :, 1) = theta(:, :, 1);
            for t = 2:n + 1
                for i = 1:n
                    for j = 1:n
                        if ~(ismember(j, SC{i}) > 0) && (theta(i, j, t) - theta(i, j, t - 1) == 1)
                            for k = 1:n
                                if ((ismember(k, SC{j}) > 0))
                                    etaz(i, k) = etax(i, k) + t-1;
                                else
                                    etaz(i, k) = etaz(i, k);
                                end
                            end
                        end
                    end


                end
                eta(:,:,t)=etaz;
                etax=etaz;
            end
        end


        %------------------
        %% Reduce the SCC accordingly

        % indexset = indices of nodes representing SCCs
        % thetaRdd = neigboring structure for the reduced graph
        %                representing SCCs
        % nuRdd = the depth of neigboring structure for the reduced graph
        %                representing SCCs

        function [indexset, thetaRdd, nuRdd] = SCCReducedStructure(G, n)

            [SC, ~] = SCCAnanyzerClass.findingSCC(G, n);

            % Convert the cell array of SCCs to a cell array of strings for comparison
            strSC = cellfun(@(x) mat2str(x), SC, 'UniformOutput', false);

            % Find unique SCCs
            uniqueStrSC = unique(strSC);

            % Convert the unique string SCCs back to cell arrays
            uniqueSCs = cellfun(@(x) eval(x), uniqueStrSC, 'UniformOutput', false);

            % Count the number of unique SCCs
            nSc = numel(uniqueSCs);

            [theta, eta] = SCCAnanyzerClass.SCCStructure(G, n);

            maxValues = zeros(1, nSc); % Initialize array to store maximum values

            for i = 1:nSc
                maxValues(i) = max(uniqueSCs{i}); % Extract maximum value from each cell array
            end
            indexset = maxValues;
            thetaRdd(:,:,1)= zeros(nSc, nSc); % Initialize the new matrix thetaRdd
            for t=1:1:n+1
                % Iterate through each unique SCC
                for i = 1:nSc
                    maxVal_i = maxValues(i);
                    iden_i=uniqueSCs{i};
                    nmi=size(iden_i,2);
                    theta_copy(:,:,:)=theta(:,:,:);

                    for iin=1:1:nmi
                        theta_copy(maxVal_i,:,t)=max(theta_copy(maxVal_i,:,t),theta_copy(iden_i(iin),:,t));

                    end
                    for iin=1:1:nmi
                        theta_copy(iden_i(iin),:,t)=theta_copy(maxVal_i,:,t);
                    end

                    for j = 1:nSc
                        % Find the maximum value of SCCs i and j

                        maxVal_j = maxValues(j);

                        % Assign the corresponding value from theta to thetaRdd
                        thetaRdd(i, j,t) = theta_copy(maxVal_i, maxVal_j,t);
                        %thetaRdd(i, j,:) = max(theta(uniqueSCs{i}, uniqueSCs{j},:));
                    end
                end
            end
            nuRdd = zeros(nSc, nSc, n + 1);
            for j = 1:nSc
                count=0;
                for i= 1:nSc
                    for k = 1:n
                        if thetaRdd(i, j, k+1) - thetaRdd(i, j, k) == 1
                            % Conditions specified in the equation
                            nuRdd(i, j, k+1) = nuRdd(i, j, k) +nnz(thetaRdd(:,j,k));
                        else
                            nuRdd(i, j, k+1) = nuRdd(i, j, k);
                        end
                    end
                end
            end
        end

    end
end


