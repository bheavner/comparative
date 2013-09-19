function [ tps fns ] = testSyntheticLethal( model )
%testSyntheticLethal test double knockouts in yeast model for synthetic
%lethal interaction
%   This function compares double knockout phenoytpe predictions with a
%   list of synthetic lethal genetic interactions downloaded from the
%   Saccharomyces Genome Database.
%
%   Input     model - a COBRA Toolbox format metabolic network
%
%   Output    tps - the number of correct "synthetic lethal" predictions - 
%               neither gene is individually essential, but the double
%               mutant can't grow.
%             fns - the number of false negative predictions - the model
%               wrongly predicts that either the double knockout can grow,
%               or one of the pair of genes is individually essential

%% acknowledgments and download instructions
% thanks to Rama Balakrishnan, of SGD, for providing these instructions for
% finding synthetic lethal interactions on SGD.
%
% "Thanks for writing in. Sure, we can build this kind of query using
% YeastMine. I have provided instructions below to retrieve all gene pairs
% of genes that interact genetically and are either 'Synthetic Rescue' or
% Synthetic lethality'.
% 
% 1) Click on the Query Builder option on the tool bar of YeastMine
% homepage (yeastmine.yeastgenome.org) 2) Click on the 'Import query from
% XML' link available in the Querybuilder box on the left corner. 3) Paste
% the following XML in the text box and hit submit.
% 
% <query name="" model="genomic" view="Interaction.gene1.primaryIdentifier
% Interaction.gene1.secondaryIdentifier Interaction.gene2.primaryIdentifier
% Interaction.gene2.secondaryIdentifier Interaction.details.annotationType
% Interaction.details.experimentType Interaction.details.phenotype
% Interaction.details.type" longDescription=""
% sortOrder="Interaction.gene1.primaryIdentifier asc">
%   <constraint path="Interaction.details.experimentType" op="="
%   value="Synthetic Lethality"/>
% </query>
% 
% 4) You will now land on a Model browser page which shows the query model.
% Click on the green Show Results page and you should get a list of
% interacting genes where experiment type = 'Synthetic Lethality'
% 
% 5) To get a list of interacting genes with Synthetic rescue type, go back
% to step 3, and paste the following XML and follow step 4.
% 
% <query name="" model="genomic" view="Interaction.gene1.primaryIdentifier
% Interaction.gene1.secondaryIdentifier Interaction.gene2.primaryIdentifier
% Interaction.gene2.secondaryIdentifier Interaction.details.annotationType
% Interaction.details.experimentType Interaction.details.phenotype
% Interaction.details.type" longDescription=""
% sortOrder="Interaction.gene1.primaryIdentifier asc">
%   <constraint path="Interaction.details.experimentType" op="="
%   value="Synthetic Rescue"/>
% </query>
% "

% After following these instructions, I processed the downloaded report
% thus:
% 
% fid=fopen('SGD synthetic lethal.txt')
% 
% A = cell(0,2);
% 
% while ~feof(fid)
%     tline = fgetl(fid);
%     tline = textscan(tline,'%s','delimiter','\t');
%     %tline is
%     %{[gene1][gene2][screentype][experiment][phenotype][interaction]}
%     A(end+1,1) = tline{1}(1); 
%     A(end,2) = tline{1}(2);
% end
% 
% save('SGD synthetic lethals','A');
% 
% fclose(fid);

    load('SGD synthetic lethals');

    KO_result = [];

    h = waitbar(0,'Checking SGD synthetic lethals ...');
    for index = 1:length(A)

        if mod(index,100) == 0
            waitbar(index/length(A),h);
        end

        if length(intersect(A(index,:),model.genes))==2
            % check each gene individually to see if they're predicted to be
            % essential
            [gene1_mutant,~,~,~] = deleteModelGenes(model,(A(index,1)));
            [gene2_mutant,~,~,~] = deleteModelGenes(model,(A(index,2)));
            [double_mutant,~,~,~] = deleteModelGenes(model,(A(index,:)));
            gene1_sln=optimizeCbModel(gene1_mutant,[],'one');
            gene2_sln=optimizeCbModel(gene2_mutant,[],'one');

            %then check the pair
            double_sln=optimizeCbModel(double_mutant,[],'one');
            KO_result(end+1,:) = ...
                [index gene1_sln.f gene2_sln.f double_sln.f];
        end
    end

    close(h);

    model_gene1_inessential = KO_result(:,2)~=0;
    model_gene2_inessential = KO_result(:,3)~=0;
    model_synthetic_lethal = KO_result(:,4)==0;

    model_inessential_genes_summed = ...
        model_gene1_inessential + model_gene2_inessential;
    model_predicts_both_genes_inessential = ...
        model_inessential_genes_summed == 2;

    % TO DOUBLE CHECK -
    % I think the correct "synthetic lethal" should be neither gene is
    % individually essential, but the double mutant can't grow. If that's
    % right, this finds the true positive rate for the model
    tps = sum(model_predicts_both_genes_inessential & model_synthetic_lethal);

    fns = length(KO_result) - tps;

end

