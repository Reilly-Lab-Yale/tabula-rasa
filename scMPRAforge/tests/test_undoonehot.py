import pandas as pd
import numpy as np
import pytest
import re

def undo_one_hot_encoding_sparse(df, reference_label="reference"):
    """Simplified sparse-friendly version that avoids idxmax."""
    

def test_undo_one_hot_encoding_sparse_simple():
    """Test the sparse-friendly undo one-hot encoding function."""
    
    # Create test data with sparse one-hot encoding
    data = {
        'C(color)[red]': [1, 0, 0, 0, 1],
        'C(color)[blue]': [0, 1, 0, 0, 0], 
        'C(color)[green]': [0, 0, 1, 0, 0],
        'C(size)[small]': [1, 0, 1, 0, 0],
        'C(size)[large]': [0, 1, 0, 0, 0],
        'other_col': [10, 20, 30, 40, 50]
    }
    
    # Convert to sparse DataFrame
    df = pd.DataFrame(data)
    for col in df.columns:
        if col.startswith('C('):
            df[col] = df[col].astype(pd.SparseDtype(int, fill_value=0))
    
    # Test the function
    result = undo_one_hot_encoding_sparse_simple(df)
    
    # Expected results
    expected_color = ['red', 'blue', 'green', 'reference', 'red']
    expected_size = ['small', 'large', 'small', 'reference', 'reference']
    
    # Assertions
    assert 'color' in result.columns
    assert 'size' in result.columns
    assert 'other_col' in result.columns
    
    # Check categorical values
    assert list(result['color']) == expected_color
    assert list(result['size']) == expected_size
    assert list(result['other_col']) == [10, 20, 30, 40, 50]
    
    # Check that one-hot columns are removed
    one_hot_cols = [col for col in df.columns if col.startswith('C(')]
    for col in one_hot_cols:
        assert col not in result.columns
    
    # Check that categorical columns are created
    assert isinstance(result['color'].dtype, pd.CategoricalDtype)
    assert isinstance(result['size'].dtype, pd.CategoricalDtype)
    
    print("✅ Basic functionality test passed")


def test_multi_hot_error():
    """Test that multi-hot encoding raises an error."""
    
    # Create invalid multi-hot data
    data = {
        'C(color)[red]': [1, 1, 0],   # Row 1 has both red and blue active
        'C(color)[blue]': [0, 1, 0],
        'other_col': [10, 20, 30]
    }
    
    df = pd.DataFrame(data)
    for col in df.columns:
        if col.startswith('C('):
            df[col] = df[col].astype(pd.SparseDtype(int, fill_value=0))
    
    # Should raise ValueError
    with pytest.raises(ValueError, match="Multi-hot encoding detected"):
        undo_one_hot_encoding_sparse_simple(df)
    
    print("✅ Multi-hot error test passed")


def test_all_reference_category():
    """Test handling of all-zero rows (reference category)."""
    
    data = {
        'C(color)[red]': [0, 0, 0],
        'C(color)[blue]': [0, 0, 0],
        'other_col': [10, 20, 30]
    }
    
    df = pd.DataFrame(data)
    for col in df.columns:
        if col.startswith('C('):
            df[col] = df[col].astype(pd.SparseDtype(int, fill_value=0))
    
    result = undo_one_hot_encoding_sparse_simple(df, reference_label="baseline")
    
    # All rows should be reference category
    assert list(result['color']) == ['baseline', 'baseline', 'baseline']
    
    print("✅ Reference category test passed")


def test_custom_reference_label():
    """Test custom reference label."""
    
    data = {
        'C(treatment)[drug_a]': [1, 0, 0],
        'C(treatment)[drug_b]': [0, 1, 0],
        'other_col': [10, 20, 30]
    }
    
    df = pd.DataFrame(data)
    for col in df.columns:
        if col.startswith('C('):
            df[col] = df[col].astype(pd.SparseDtype(int, fill_value=0))
    
    result = undo_one_hot_encoding_sparse_simple(df, reference_label="placebo")
    
    expected = ['drug_a', 'drug_b', 'placebo']
    assert list(result['treatment']) == expected
    
    print("✅ Custom reference label test passed")


if __name__ == "__main__":
    # Run tests
    test_undo_one_hot_encoding_sparse_simple()
    test_multi_hot_error() 
    test_all_reference_category()
    test_custom_reference_label()
    
    print("\n🎉 All tests passed!")