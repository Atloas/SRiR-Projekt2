clc;
clf;
clear;

dataString = fileread('resultdata.txt');
fullData = sscanf(dataString(10:end), '%d;%f;%f;%f\n', [4 Inf]);

trackedBodies = max(fullData(1, :)) + 1;
datapointCount = size(fullData, 2)/trackedBodies;

fullData = reshape(fullData, 4, 36, datapointCount);

myVideo = VideoWriter('visualization');
myVideo.FrameRate = 30;
open(myVideo)
observedBody = 15;
windowSize = 4e9;

x0 = 10;
y0 = 10;
width = 505;
height = 480;

sun = 1;
planets = [2 3 4 6 10 15 24 30];
moonsAndAsteroids = [5 7 8 9 11 12 13 14 16 17 18 19 20 21 22 23 25 26 27 28 29 31 32 33 34];
probes = [35 36];

for i = 1:9:datapointCount
    scatter3(fullData(2, probes, i), fullData(3, probes, i), fullData(4, probes, i), 'filled', 'k');
    hold on;
    scatter3(fullData(2, moonsAndAsteroids, i), fullData(3, moonsAndAsteroids, i), fullData(4, moonsAndAsteroids, i), 'filled', 'b');
    scatter3(fullData(2, planets, i), fullData(3, planets, i), fullData(4, planets, i), 'filled', 'g');
    scatter3(fullData(2, sun, i), fullData(3, sun, i), fullData(4, sun, i), 'filled', 'y');
    hold off;
    
    xlim([fullData(2, observedBody, i) - windowSize, fullData(2, observedBody, i) + windowSize]);
    ylim([fullData(3, observedBody, i) - windowSize, fullData(3, observedBody, i) + windowSize]);
    zlim([fullData(4, observedBody, i) - windowSize, fullData(4, observedBody, i) + windowSize]);
    view([0 0 1]);
    title('Uk³ad planetarny Saturna');
    xlabel('X [m]');
    ylabel('Y [m]');
    
    pause(0.005);
    set(gcf, 'position', [x0, y0, width, height]);
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo)
