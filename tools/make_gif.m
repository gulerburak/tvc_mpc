function make_gif(avi_file, gif_file, max_mb)
    % MAKE_GIF  Convert the simulation AVI into a README-friendly animated GIF.
    %
    %   make_gif()
    %   make_gif(avi_file, gif_file, max_mb)
    %
    % GitHub will not play an AVI inline in a README, so the animation written by
    % run_mpc_sim.m is transcoded to a looping GIF here. Frames are downscaled and
    % decimated, and a single shared 128-colour palette is derived from samples
    % across the whole clip so the colours stay stable frame to frame.
    %
    % Successively more aggressive settings are tried until the result fits in
    % max_mb. This is an authoring step, not part of the simulation: it is the one
    % file in this repo that needs the Image Processing Toolbox (rgb2ind,
    % imresize), and you only need it if you are regenerating the asset.

    root = fileparts(fileparts(mfilename('fullpath')));

    if nargin < 1 || isempty(avi_file); avi_file = fullfile(root, 'tvc_animation.avi'); end
    if nargin < 2 || isempty(gif_file); gif_file = fullfile(root, 'assets', 'tvc_animation.gif'); end
    if nargin < 3 || isempty(max_mb); max_mb = 5; end

    assert(isfile(avi_file), 'AVI not found: %s (run main_tvc_mpc first)', avi_file);

    asset_dir = fileparts(gif_file);
    if ~exist(asset_dir, 'dir'); mkdir(asset_dir); end

    % width [px], frame stride — tried in order until the file fits
    settings = [900 2; 800 3; 720 3; 640 4];

    for s = 1:size(settings, 1)
        width = settings(s, 1);
        stride = settings(s, 2);

        write_gif(avi_file, gif_file, width, stride);

        d = dir(gif_file);
        mb = d.bytes / 1048576;
        fprintf('  %4d px / every %d frames -> %.2f MB\n', width, stride, mb);

        if mb <= max_mb
            fprintf('GIF written: %s  (%.2f MB)\n', gif_file, mb);
            return;
        end

    end

    warning('make_gif:tooBig', 'Could not get under %g MB; kept the smallest attempt.', max_mb);
end

function write_gif(avi_file, gif_file, width, stride)
    vr = VideoReader(avi_file);
    idx = 1:stride:vr.NumFrames;
    n = numel(idx);

    scale = width / vr.Width;
    h = round(vr.Height * scale);

    % Pass 1: resize, keeping a subsample to derive the shared palette from.
    small = zeros(h, width, 3, n, 'uint8');

    for i = 1:n
        small(:, :, :, i) = imresize(read(vr, idx(i)), [h width]);
    end

    n_probe = min(n, 10);
    probe = small(:, :, :, round(linspace(1, n, n_probe)));
    probe = reshape(permute(probe, [1 4 2 3]), h * n_probe, width, 3);
    [~, map] = rgb2ind(probe, 128, 'nodither');

    % Pass 2: index every frame against that one palette.
    frames = zeros(h, width, 1, n, 'uint8');

    for i = 1:n
        frames(:, :, 1, i) = rgb2ind(small(:, :, :, i), map, 'nodither');
    end

    delay = stride / vr.FrameRate;
    imwrite(frames, map, gif_file, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
end
