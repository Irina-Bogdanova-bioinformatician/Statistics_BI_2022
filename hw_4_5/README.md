# Differential gene expression assessment

**homework_lecture_5.ipynb**
Colab notebook with some practice on differential gene expression assessment.

**compare_expressions.py**

Compare genes expressions between two cell types.\
Input:\
--first_cell_type_expressions_path: Path to the table (csv file) with gene expressions of first cell type (required!)\
--second_cell_type_expressions_path: Path to the table (csv file) with gene expressions of second cell type (required!)\
--save_results_table: Path to the output table with gene expressions comparison results (default="expression_comparison_results.csv")

Output:\
You get a table (.csv file) with confidence intervals test results, z-test results, z-test p-values, and a difference between mean expressions.

### Requirements
python>=3.6

### Create and activate virtual environment
In work directory run:
~~~sh
python3 -m venv new_venv_name
source new_venv_name/bin/activate
~~~

### Install requirements from requirements.txt
Run:
~~~sh
pip install -r requirements.txt
~~~

### Run compare_expressions.py
Run:
~~~sh
python3 compare_expressions.py \
--first_cell_type_expressions_path files_for_test/b_cells_expression_data.csv \
--second_cell_type_expressions_path files_for_test/nk_cells_expression_data.csv \
--save_results_table path_to_your_awesome_results.csv
~~~
