import pytest
from . import create_input_template_for_scoring_module

def test__scoring__baseline_macc_score_run():
    df = create_input_template_for_scoring_module()
    df.to_csv('/Users/cranneyc/Desktop/df.csv', index=False)