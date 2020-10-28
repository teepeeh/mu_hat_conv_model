



%% Plot just stim & Response
IR = IR_wrapper.IR_designmatrix;  


x = cell2mat(IR.stimuli.stimuli);
x2 = cell2mat(IR.response.response);

close all
a = plot(x(:,1),'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
b = plot(x(:,2),'.', 'LineWidth',2,...
    'MarkerSize',10);
c = plot(x2(:,1)+0.1,'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
d = plot(x2(:,2)+0.1,'.', 'LineWidth',2,...
    'MarkerSize',10);

ylim([0.98 1.12])

legend([a b c d],{'Blue Left', 'Orange Left','Resp L', 'Resp R'})
set(gcf,'color','white')


%% Plot OUTCOME stim & Response
IR = IR_wrapper.IR_designmatrix;  


x = cell2mat(IR.stimuli.stimuli);
x2 = cell2mat(IR.response.response);
x3 = cell2mat(IR.outcomes.outcomes);

close all
a = plot(x(:,1),'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
b = plot(x(:,2),'.', 'LineWidth',2,...
    'MarkerSize',10)
c = plot(x2(:,1)+0.02,'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
d = plot(x2(:,2)+0.02,'.', 'LineWidth',2,...
    'MarkerSize',10);

e = plot(x2(:,1)+0.04,'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
f = plot(x2(:,2)+0.04,'.', 'LineWidth',2,...
    'MarkerSize',10);

ylim([0.99 1.05])

legend([a b c d e f],{'Blue Left', 'Orange Left','Resp L', 'Resp R','Win','Lose'})
set(gcf,'color','white')





%% Plot OUTCOME ONLY

IR = IR_wrapper.IR_designmatrix;  


x = cell2mat(IR.outcomes(1));


close all
a = plot(x(:,1),'.', 'LineWidth',2,...
    'MarkerSize',10);
hold on
b = plot(x(:,2),'.', 'LineWidth',2,...
    'MarkerSize',10)

ylim([0.99 1.05])

legend([a b],{'Win','Lose'})
set(gcf,'color','white')

