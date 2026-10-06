PYTHON ?= python3
PORT ?= 8501

.PHONY: setup pipeline dashboard clean

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_pipeline.py

dashboard:
	$(PYTHON) -m streamlit run dashboard/app.py \
		--server.address 0.0.0.0 --server.port $(PORT) \
		--server.headless true --browser.gatherUsageStats false

clean:
	rm -rf cell_counts.db cell_counts.db.tmp outputs
