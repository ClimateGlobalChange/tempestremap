// This is free and unencumbered software released into the public domain.
// 
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
// 
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// 
// For more information, please refer to <http://unlicense.org>
//
// Code obtained from:
//   Ebbeke, Tim (2021) "interval-tree" https://github.com/5cript/interval-tree,
//   accessed August 22, 2021 commit f5bff66f88daf47b03a5a612c58f09e2a399c199.
//

#pragma once

namespace lib_interval_tree
{
    template <typename numerical_type, typename interval_kind_, typename data_type_>
    struct interval;

    template <typename IntervalT>
    class interval_tree;

    template <typename numerical_type, typename data_type, typename interval_type>
    class node;

    template <typename node_type, typename owner_type>
    class basic_interval_tree_iterator;

    template <typename node_type>
    class const_interval_tree_iterator;

    template <typename node_type>
    class interval_tree_iterator;
}
