from config_code.read import main as run_data_checks


def main():
    print("HPV16 and Human lncRNA Complementarity Project")
    print("Starting data quality checks...\n")

    run_data_checks()

    print("\nData quality checks complete.")


if __name__ == "__main__":
    main()
