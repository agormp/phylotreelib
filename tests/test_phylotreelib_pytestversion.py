import phylotreelib as pt
import pytest

###################################################################################################
###################################################################################################

# Tests for loose functions

###################################################################################################
###################################################################################################

class Test_remove_comments:

    def test_balanced(self, balanced_input):
        assert pt.remove_comments(balanced_input) == 'Hello  World!\n'

    def test_unbalanced(self, unbalanced_input):
        with pytest.raises(pt.TreeError):
            pt.remove_comments(unbalanced_input)

    def test_nested_comments(self, nested_comments):
        assert pt.remove_comments(nested_comments) == 'Hello  World!'