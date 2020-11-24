%% Simulated RF mapping experiment.
% A linear, single-opponent neuron with an L-M Gaussian RF is probed with
% simultaneous cone-isolating contrast increments and decrements.
% A faint cone-opponent halo appears outside of the RF because the STA
% must contain equal parts contrast increments and decrements.

% Simulated L-cone map
sigma = .5; % Standard deviation of single-opponent (L-M) RF
bounds = 3;
boxwidth = .5;
niter = 2000;

data = [];
for i = 1:niter
    inc_pos = unifrnd(-bounds,bounds,2,1);
    dec_pos = unifrnd(-bounds,bounds,2,1);
    while abs(inc_pos(1)-dec_pos(1))<boxwidth & abs(inc_pos(2)-dec_pos(2))<boxwidth
        dec_pos = unifrnd(-bounds,bounds,2,1); % No overlapping increments with decrements
    end
    r_inc_x = normcdf(inc_pos(1)+[-.5 .5]*boxwidth, 0, sigma);
    r_inc_y = normcdf(inc_pos(2)+[-.5 .5]*boxwidth, 0, sigma);
    r_inc = diff(r_inc_x)*diff(r_inc_y);
    r_dec_x = normcdf(dec_pos(1)+[-.5 .5]*boxwidth, 0, sigma);
    r_dec_y = normcdf(dec_pos(2)+[-.5 .5]*boxwidth, 0, sigma);
    r_dec = diff(r_dec_x)*diff(r_dec_y);
    r = r_inc-r_dec;
    data = [data; r inc_pos' dec_pos'];
end

% Calculating STA
[x,y] = meshgrid(linspace(-bounds, bounds,20));
inc_map = zeros(size(x));
dec_map = zeros(size(x));
for i = 1:numel(x)
    L_inc = abs(data(:,2)-x(i)) < boxwidth & abs(data(:,3)-y(i)) < boxwidth;
    L_dec = abs(data(:,4)-x(i)) < boxwidth & abs(data(:,5)-y(i)) < boxwidth;
    inc_map(i) = (data(:,1)'*L_inc);
    dec_map(i) = (data(:,1)'*L_dec);
end

% For image rendering, putting min at 0 and max at 255
min_val = min([inc_map(:); dec_map(:)]);
max_val = max([inc_map(:); dec_map(:)]);

% actual RF, normalized to [0:max_val]
RF = normpdf(x, 0, sigma).*normpdf(y, 0, sigma);
RF = RF./max(RF(:))*max_val;

figure; 
colormap(jet(255))
subplot(2,2,1);
image(255*(RF-min_val)./(max_val-min_val)); 
axis square; set(gca,'Xtick',[],'Ytick',[]);
subplot(2,2,2);
image(255*(-RF-min_val)./(max_val-min_val)); 
axis square; set(gca,'Xtick',[],'Ytick',[]);

subplot(2,2,3);
image(255*(inc_map-min_val)./(max_val-min_val)); 
axis square; set(gca,'Xtick',[],'Ytick',[]);
subplot(2,2,4);
image(255*(dec_map-min_val)./(max_val-min_val)); 
axis square; set(gca,'Xtick',[],'Ytick',[]);
figure; axes; hold on;
plot(x(1,:),inc_map(10,:),'r-');
plot(x(1,:),dec_map(10,:),'b-');
