rule test:
    output:
        "test.txt"
    shell:
        """
        touch test.txt
        """
