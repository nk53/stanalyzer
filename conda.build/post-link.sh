if [ "$PKG_NAME" != "stanalyzer-dev" ]; then
    exit 0
fi

cat << END >> $PREFIX/.messages.txt


WARNING: DEVELOPMENT INSTALLATION IS NOT COMPLETE.

Complete stanalyzer installation with the following command:
    sta-finish-install

END
