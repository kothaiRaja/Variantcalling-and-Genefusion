process CHECK_JAVA {
    tag "Check Java"
    container null

    output:
    path "java_check.txt"

    script:
    """
    # Check if Java is installed
    if ! java -version &>/dev/null; then
        echo "Java is not installed or not in PATH." > java_check.txt
        exit 1
    fi

    # Extract the Java version
    JAVA_VERSION=\$(java -version 2>&1 | awk -F '"' '/version/ {print \$2}' | awk -F '.' '{print \$1}')

    # Check if Java version is 21 or higher
    if [ "\$JAVA_VERSION" -lt 21 ]; then
        echo "Java version \$JAVA_VERSION detected. Update to Java 21 or higher is required." > java_check.txt
        exit 1
    else
        echo "Java version \$JAVA_VERSION is available and meets the requirements." > java_check.txt
    fi
    """
}
