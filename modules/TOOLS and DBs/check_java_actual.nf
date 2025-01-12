process CHECK_JAVA {
    tag "Check Java"
	container null
    output:
    path "java_check.txt"

    script:
    """
    if ! java -version &>/dev/null; then
        echo "Java is not installed or not in PATH." > java_check.txt
        exit 1
    else
        echo "Java is available." > java_check.txt
    fi
    """
}