function nodes = collect_subtree(parent, root)
    nodes = root;
    children = find(parent == root);
    nodes = [nodes, children];
    % for c = children(:).'
    %     nodes = [nodes, collect_subtree(parent, c)]; %#ok<AGROW>
    % end
end