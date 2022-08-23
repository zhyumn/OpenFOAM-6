#include<vector>
#include<fstream>
struct temp_node
{
    int left, right;
    void* ptr;
};
template <class T>
void traverse
(
    T* node_ptr, int node_label,
    std::vector<temp_node>& node, std::vector<void*>& leaf,
    int& node_tail, int& leaf_tail
)
{
    node[node_label].ptr = node_ptr;
    if (left_is_node)
    {
        node[node_label].left = node_tail++;
        traverse(node_ptr->left_node, node[node_label].left, node, leaf, node_tail, leaf_tail);
    }
    else
        (
            leaf[leaf_tail] = node_ptr->left_leaf;
    node[node_label].left = -(++leaf_tail);
    )

        if (right_is_node)
        {
            node[node_label].right = node_tail++;
            traverse(node_ptr->right_node, node[node_label].right, node, leaf, node_tail, leaf_tail);
        }
        else
            (
                leaf[leaf_tail] = node_ptr->right_leaf;
    node[node_label].right = -(++leaf_tail);
    )

}

template <class T>
void treewriter(T tree)
{
    int leafsize = tree.size();
    if (leafsize <= 0)
        return;
    int nodesize = leafsize - 1;
    int node_tail = 1;
    int leaf_tail = 0;
    std::vector<temp_node> node_(nodesize);
    std::vector<void*> leaf_(leafsize);
    if (leafsize == 1)
    {
        leaf_[0] = tree.root().leftleaf();
    }
    traverse(tree.root(), 0, node_, leaf_, node_tail, leaf_tail);

    std::ofstream nodewriter("node.txt");
    std::ofstream leafwriter("leaf.txt");
    for (int i = 0;i < nodesize;i++)
    {
        nodewriter << node_[i].left << " " << node_[i].right << " ";
        node_[i].ptr->write(nodewriter);
        nodewriter << endl;
    }
    for (int i = 0;i < leafsize;i++)
    {
        leaf_[i].ptr->write(leafwriter);
        leafwriter << endl;
    }

}