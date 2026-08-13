function save_figure(fig, name)
    % SAVE_FIGURE  Export a figure twice, once per audience.
    %
    %   save_figure(fig, name)
    %
    %   fig   figure handle (pass gcf for the current figure)
    %   name  base filename, WITHOUT extension
    %
    % Writes a vector PDF to report/figures/ for the write-up, and a 150 dpi
    % PNG to assets/ for the README. Both directories are created on demand and
    % are resolved relative to this file, so the output lands next to the code
    % regardless of the current working directory.

    root = fileparts(mfilename('fullpath'));
    fig_dir = fullfile(root, 'report', 'figures');
    asset_dir = fullfile(root, 'assets');

    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
    if ~exist(asset_dir, 'dir'); mkdir(asset_dir); end

    exportgraphics(fig, fullfile(fig_dir, [name '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, fullfile(asset_dir, [name '.png']), ...
        'Resolution', 150, 'BackgroundColor', 'white');
end
