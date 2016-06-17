function [dataset, Y_train, Y_test, label_train, label_test] = train_test_split(...
    dataset, N_train)
myrng();
fprintf('Picking Train and Test set ...');
switch dataset
        case 'myARgender'
            AR_gender_fn = fullfile('data', 'myARgender.mat');
            load(AR_gender_fn);            
            Y_train = normc(double(Y_train));
            train_range = label_to_range(label_train);
            Y_train = PickDfromY(Y_train, train_range, N_train);
            C = numel(train_range) - 1;
            label_train = range_to_label(N_train*(0:C));
            Y_test = normc(double(Y_test));

        case 'myARreduce'
            dataset = 'test mode';
            load('data/AR_EigenFace.mat');
            % ---------------  -------------------------
            Y_train = normc(tr_dat);
            Y_test = normc(tt_dat);
            % ---------------  -------------------------
            label_train = trls;
            label_test = ttls;
            
        case 'myFlower'
            dataset = 'myFlower102';
            % load(fullfile('data', dataset);
            load(fullfile('data',strcat(dataset, '.mat')));

            Y_train = normc(double(Y_train));
            train_range = label_to_range(label_train);
            Y_train = PickDfromY(Y_train, train_range, N_train);
            C = numel(train_range) - 1;

            label_train = range_to_label(N_train*(0:C));

            Y_test = normc(double(Y_test));
        otherwise
            [Y_train, label_train, Y_test, label_test] = ...
                pickTrainTest_2(dataset, N_train);
    end
    fprintf('DONE\n');
end 
