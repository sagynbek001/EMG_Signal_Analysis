% This script analyzes the correlation between the Force and EMG signals
%{
Name: Sagynbek Talgatuly, Student Number: st4121			 	
Date: December 16. 2020							 
Program: Assignment5.m							 	 
Description: 
    This program works with two data signals: EMG signal and
    Force signal. It filters them, crops them, normalizes them.
    The program, then, identifies equal segments of the data 
    and calculates the average segment. Lastly, it plots the 
    relation graph of this two averages and finds the best fit 
    function for the period of contraction.
%}
%% Task 0
% clearing the workspace and command window
clc;
clear;
%% Task 1

% loading the data file
load('Data.mat');

% plotting the graphs of the original data
figure(1);
subplot (2,1,1);
plot (timestamps, EMG_signal);
title ('EMG signal (original)');
xlabel ('Time (seconds)');
ylabel ('EMG (in volts)');
subplot (2,1,2);
plot (timestamps, Force_signal);
title ('Force signal (original)');
xlabel ('Time (seconds)');
ylabel ('Force');

% applying the low-pass filter with a cutoff frequency 200Hz
% sampling frequency is given in the problem statement
% cutoff frequency is given in the problem statement
sampling_frequency = 1000;
cutoff_frequency = 200;
EMG_signal = lowpass(EMG_signal,cutoff_frequency,sampling_frequency);
Force_signal = lowpass(Force_signal,cutoff_frequency,sampling_frequency);

% plotting the graphs of the filtered data
figure(2);
subplot (2,1,1);
plot (timestamps, EMG_signal);
title ('EMG signal (filtered)');
xlabel ('Time (seconds)');
ylabel ('EMG (in volts)');
subplot (2,1,2);
plot (timestamps, Force_signal);
title ('Force signal (filtered)');
xlabel ('Time (seconds)');
ylabel ('Force');
clear sampling_frequency;
clear cutoff_frequency;

% cropping the endings of the data
length_data = length(timestamps);

EMG_signal = EMG_signal(501:length_data-250);
Force_signal = Force_signal(501:length_data-250);
timestamps = timestamps(501:length_data-250);

length_data = length(timestamps);

% normalizing the force signal
global_max_of_f = Force_signal(1);
global_min_of_f = Force_signal(1);

% finding the global max and min
for i = 1:length_data
    if Force_signal(i) > global_max_of_f
        global_max_of_f = Force_signal(i);
    end
    if Force_signal(i) < global_min_of_f
        global_min_of_f = Force_signal(i);
    end
end

% applying the normalization formula
for j = 1:length_data
    top = (Force_signal(j) - global_min_of_f);
    bot = (global_max_of_f - global_min_of_f);
    Force_signal(j) = 100 * top / bot;
end

clear top;
clear bot;
clear i;
clear j;
clear global_max_of_f;
clear global_min_of_f;

% plotting the graphs of cropped and normalized data
figure(3);
subplot (2,1,1);
plot (timestamps, EMG_signal);
title ('EMG signal (cropped)');
xlabel ('Time (seconds)');
ylabel ('EMG (in volts)');
subplot (2,1,2);
plot (timestamps, Force_signal);
title ('Force signal (cropped and normalized)');
xlabel ('Time (seconds)');
ylabel ('Force (in percentage of MVC)');
%% Task 2 & 3

% applying the rms envelope function for EMG signal
window_length = 150;
[upper_envelope, bottom_envelope] = envelope((EMG_signal),window_length,'rms');
clear window_length;

% setting the length of one cycle
% creating the matrix for the average segment
segment_length = 500;
average_segment = zeros(2,segment_length);

maxPoints_Force = zeros(31);
maxPoints_EMG = zeros(31);

% finding the peaks of the Force signal
i = 1;
for k = 1:31
    while 1==1
        if (Force_signal(i) > 25 && Force_signal(i) < Force_signal(i+1))
            % if some point is greater than 25 and at this point the graph
            % is increasing, then we break the loop at this point
            break;
        end
        i = i + 1;
    end
    
    MAX = i + 1;
    
    j = 1;
    % from that point it starts looking for the peak 
    % it seeks for the peak until the point is less than 25
    while Force_signal(i+j) > 25
        if Force_signal(i+j) > Force_signal(MAX)
            MAX = i+j;
        end
        j = j + 1;
    end
    
    i = i + j;
    % the peak is found and assigned to the vector
    % this process goes 31 times
    maxPoints_Force(k) = MAX;
end
clear i;
clear j;
clear MAX;

% finding the peaks of the EMG signal
k = 1;
for i = 250:length_data-350
    % in this loop, all the points from 250 to 28900 are analyzed
    % for this algorith, the program uses a flag variable
    flag = 1;
    % if the current element is less than 0.19 then flag is set to 0
    % because we know that the lowest peak is greater than 0.2 (manually checked)
    if upper_envelope(i)<0.19
        flag = 0;
    end
    % if flag is still 1, the next condition is checked
    % the condition that this point is greater than all the points in a
    % certain scope
    if flag == 1
        for j = 1:249
            if upper_envelope(i)<upper_envelope(i+j) || upper_envelope(i)<upper_envelope(i-j)
                flag = 0;
            end
        end
    end
    % if the point passed the two conditions above
    % then it concludes that this point is indeed the peak
    if flag == 1
        maxPoints_EMG(k) = i;
        k = k + 1;
    end
end

clear flag;
clear i;
clear j;

% calculating the average segment
for k = 1:31
    average_segment(1,:) = average_segment(1,:) + Force_signal(maxPoints_Force(k)-(segment_length/2-1):maxPoints_Force(k)+segment_length/2)';
    average_segment(2,:) = average_segment(2,:) + upper_envelope(maxPoints_EMG(k)-(segment_length/2-1):maxPoints_EMG(k)+(segment_length/2))';
end
clear k;
clear maxPoints;
average_segment = average_segment / 31;

% plotting the average segment graphs of EMG and Force signals
figure(4);
subplot(3,1,1);
plot(1:segment_length,average_segment(2,:));
title ('Average segment for EMG');
xlabel ('Time (milliseconds)');
ylabel ('EMG (in volts)');
subplot(3,1,2);
plot(1:segment_length,average_segment(1,:));
title ('Average segment for Force');
xlabel ('Time (milliseconds)');
ylabel ('Force (MVC)');
subplot(3,1,3);
plot(average_segment(2,(1:segment_length/2)),average_segment(1,(1:segment_length/2)));
hold on;
plot(average_segment(2,(segment_length/2:segment_length)),average_segment(1,(segment_length/2:segment_length)));
%% Task 4

% finding the best fit function
polynom_f_coefficients = polyfit(average_segment(2,(1:segment_length/2)),average_segment(1,(1:segment_length/2)),1);
best_fit_function = polyval(polynom_f_coefficients,average_segment(2,(1:segment_length/2)));
% plotting the best fit function on the same figure 
plot(average_segment(2,(1:segment_length/2)), best_fit_function);
title ('Correlation graph');
xlabel ('EMG (in volts)');
ylabel ('Force (MVC)');
legend({'Contraction', 'Release', 'Best fit line'}, 'location', 'northwest');
hold off;