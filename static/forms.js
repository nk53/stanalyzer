/**
 * Given a string that looks like ``some[text][with][brackets]``, return the
 * part inside the first set of square brackets `[]`
 *
 * @param {String} str - string to check
 * @param {Number} n - index of key to return
 * @returns {String} key
 */
function nth_key(str, n=0) {
    const subscript_re = /\[([^\]]*)\]/g; // match text inside these: []

    // each match looks like ["[text]", "text"], so matches[0][1] gets the
    // first subscript, not including brackets
    const matches = [...str.matchAll(subscript_re)];
    const key = matches[n][1];

    return key;
}

/**
 * Shortcut for :js:func:`nth_key` with ``n=0``.
 */
function first_key(str) {
    return nth_key(str);
}

/**
 * Shortcuts to select jQuery elements by name or other attribute
 */

/**
 * Return jQuery object with name selector; for name-only selections
 *
 * @param {String} sel - name to select
 * @param {String} rel - relational operator
 * @returns {jQuery} elements
 */
function jq_name(sel, rel='=') {
    return $(sel_name(sel, rel));
}

/**
 * Return an name selector string; for use in compound selections
 *
 * @returns {String} CSS selector
 */
function sel_name(sel, rel='=') {
    return `#${active_page} [name${rel}"${sel}"]`;
}

/**
 * Select elements by a single attribute
 *
 * @returns {jQuery} elements
 */
function jq_attr(attr, sel, rel='=') {
    return $(sel_attr(attr, sel, rel));
}

/**
 * Return an attribute selector string
 *
 * @returns {String} CSS selector
 */
function sel_attr(attr, sel, rel='=') {
    return `#${active_page} [${attr}${rel}"${sel}"]`;
}

/**
 * Return an ID selector that is unique across all pages
 *
 * @returns {String} CSS selector
 */
function sel_id(id, rel='=') {
    if (['=', '^='].includes(rel))
        return `[id${rel}"${active_page}_${id}"]`;
    return `[id${rel}"${id}"]`;
}

/**
 * Return a jQuery objected selected by ID (unique across pages)
 *
 * @returns {jQuery} element
 */
function jq_id(id, rel='=') {
    return $(sel_id(id, rel));
}

/**
 * Select an input's label
 *
 * @param {String} for_id - ID of ``<input>`` element
 * @param {String} [rel]
 * @returns {jQuery} element
 */
function jq_label(for_id, rel='=') {
    return $(sel_label(for_id, rel));
}

/**
 * @param {String} for_id - ID of ``<input>`` element
 * @param {String} [rel]
 * @returns {String} CSS selector
 */
function sel_label(for_id, rel='=') {
    return `label[for${rel}"${active_page}_${for_id}"]`;
}

/**
 *  Callback functions should take 2 args:
 *      target  the element to update
 *      src     the form input whose value is used for updating
 */

/**
 * Converts a number to a comma separated list of numbers from 1 to the
 * target's number. E.g., ``4`` becomes ``1,2,3,4``.
 */
function num_to_range_list(target, src) {
    let from = $(target).attr('data-range-start');
    if (!from)
        from = 0;

    const len = Number.parseInt( src.val() ) - from;

    let nums = Array(len);
    for (let i = 0; i < len; ++i)
        nums[i] = from++;

    target.html(nums.join(','));
}

/**
 * Applies updates from shell script template to the code textarea
 */
function update_code() {
    jq_name('code').val( jq_id('script_tpl').text() );
}

/**
 * Returns a callback function to toggle visibility of an email option
 * subsection.
 *
 * @params {String} assign_page - page_id of the currently active page
 * @return {function} callback - the toggler
 */
function toggle_email(assign_page) {
    var assigned_page = assign_page;
    if (!assign_page)
        console.log('Missing page assignment for toggle_email()');

    let box_sel = `#${assigned_page} [name='PBS[mail][use]']`;
    var box = $(box_sel);
    var email_opts = $(`#${assigned_page} .email-options`);

    if (!box.length)
        console.log(`Missing box: ${box_sel}`);

    box = box[0];

    return function() {
        if (box.checked)
            email_opts.show();
        else
            email_opts.hide();
    }
}

function update_visibility(elem) {
    function update_elem_visibility(elem) {
        const elem_id = elem.attr('id');
        const elem_label = $(`label[for="${elem_id}"]`);
        const elem_type = elem.attr('type');

        let is_allowed = true;

        const new_rules = elem.data('visibility');
        for (const [rule_str, rule_value] of Object.entries(new_rules)) {
            const rule = rule_str.split('_');
            const rule_type = rule[0];

            const ref_elem_sel = rule.slice(1).join('_');
            const ref_elem = $(`#${active_page} [name="${ref_elem_sel}"]`);
            const ref_elem_type = ref_elem.attr('type');

            let vis_init = ref_elem.data('vis-init');
            if (!vis_init)
                vis_init = {};
            if (!vis_init[elem_id]) {
                ref_elem.on('change', () => update_elem_visibility(elem));
                vis_init[elem_id] = 1;
            }
            ref_elem.data('vis-init', vis_init);

            let ref_value;
            switch (ref_elem_type) {
                case 'radio':
                    ref_value = ref_elem.filter(':checked').val();
                    break;
                default:
                    ref_value = ref_elem.val();
            }

            switch (rule_type) {
                case 'allowed':
                    is_allowed = is_allowed && (ref_value == rule_value);
                    break;
                case 'disallowed':
                    is_allowed = is_allowed && (ref_value != rule_value);
                    break;
                default:
                    console.log(`Unrecognized rule type: '${rule_type}'`);
            }

            if (!is_allowed)
                break;
        }

        let target = elem;
        if (elem_type == 'button')
            // hide button's div.ui-field-contain parent
            target = $(target.parents()[1]);

        if (is_allowed) {
            target.show();
            elem_label.show();
            target.filter('textarea').trigger('keyup');
        } else {
            target.hide();
            elem_label.hide();
        }
    }

    if (elem)
        update_elem_visibility($(elem));
    else
        $(`#${active_page} [data-visibility]`).each(function(index, elem) {
          update_elem_visibility($(elem));
        });
}

function toggle_visibility(elem) {
    function toggle_elem_visibility(toggle) {
        toggle = $(toggle);
        const target = toggle.data('target');
        if (toggle[0].checked)
            target.show();
        else
            target.hide();
    };

    if (elem)
        toggle_elem_visibility(elem);
    else
        $(`#${active_page} [data-target]`).each(
            (idx, elem) => toggle_elem_visibility(elem)
        );
}

function get_current_project(elem=null) {
    if (elem === null)
        elem = jq_id('project');

    const selected = Number.parseInt($(elem).val());

    if (!Number.isInteger(selected))
        return undefined;

    const settings = JSON.parse(jq_id('json').val());
    const selected_project = settings[selected];

    return selected_project;
}

function confirm_project_then(elem, callback) {
    const message = $(elem).data('confirm-msg');

    if (confirm(message))
        return callback(elem);

    return false;
}

function submit_async(elem=null) {
    const form = $(`form#${active_page}_form`);

    const settings = {
        url: form.attr('action'),
        method: form.attr('method'),
        data: JSON.stringify(form.serializeJSON()),
        contentType: 'application/json',
        success: console.log,
        failure: console.log,
    }

    const eval_settings = ["data", "success", "failure"];

    if (elem !== null) {
        elem = $(elem);
        for (let key of Object.keys(settings)) {
            const value = elem.data(key);
            if (value !== undefined) {
                if (eval_settings.includes(key))
                    settings[key] = eval(value);
                else
                    settings[key] = value;
            }
        }
    }

    console.log(settings);
    $.ajax(settings);
}

function update_projects_menu(response_json) {
    const request_type = response_json.request_type;
    const new_menu = response_json.data;
    const settings_elem = jq_id('json');
    const menu_elem = jq_name('project');

    // update in-form JSON representation
    settings_elem.val(JSON.stringify(new_menu));

    // TODO: on delete or post, set active menu item
    // what to do if only remaining item is add_new?
    //if (request_type == 'post') {

    //}

    // remove previous options
    menu_elem.children(':not(option[value=add_new])').remove();
    $.each(new_menu, function(id, obj) {
        const text = obj.title;
        menu_elem.prepend(`<option value="${id}">${text}</option>`);
    });
}

function load_project_json(elem) {
    const selected = Number.parseInt($(elem).val());
    if (Number.isInteger(selected)) {
        const settings = JSON.parse(jq_id('json').val());
        const selected_project = settings[selected];

        if (!selected_project)
            console.log("No selected project");

        for (const [key, value] of Object.entries(selected_project)) {
            if (key == 'time_step') {
                const [ts_num, ts_scale] = value.split(' ');
                jq_name('time_step[num]').val(ts_num);
                jq_label(`time_step[scale]-${ts_scale}`).click();

                continue;
            }
            const update_elem = jq_name(key);
            if (update_elem) {
                if (update_elem.attr('type') == 'radio')
                    jq_label(`${key}-${value}`).click();
                else
                    update_elem.val(value);
            }
        }
    } else {
        $(`#${active_page} form :is(input[type=text], input[type=number], textarea)`).each(
            function (index, elem) {
                elem = $(elem);
                const name = elem.attr('name');
                console.log(name);
                if (name && name.includes('time_step'))
                    console.log('hi');
                elem.val(elem.data('default') || '');
                elem.filter('textarea').trigger('keyup');
            });
    }
}

var active_page = null;
var loaded_pages = [];

/**
 * Handles page changes, initializes event handlers
 *
 * @param {jQuery.Event} event - this is ignored
 * @param {Object} ui - see jQuery mobile `pagecontainerload`_
 */
function setup_form(event, ui) {
    active_page = ui.toPage[0].id;

    if (loaded_pages.includes(active_page))
        return;

    if (['home', '404'].includes(active_page))
        return;

    $(`#${active_page} [data-toggle]`).each(function(index, elem) {
        elem = $(elem);
        const toggle_id = elem.attr('data-toggle');
        const toggle = $(`#${active_page} #${toggle_id}`);

        toggle.data('target', elem);
        toggle.on('click', () => toggle_visibility(toggle));
        toggle_visibility(toggle);
    });

    $(`#${active_page} [data-setup]`).each(function(index, elem)  {
        eval($(elem).data('setup'));
    });

    toggle_visibility();
    update_visibility();
}

$(document).on('pageload', function(event, ui) {
    active_page = ui.toPage[0].id;
    let index = loaded_pages.indexOf(active_page);
    if (index !== -1)
        loaded_pages.splice(index, 1); // remove it to allow re-loading
});
$(document).on('pagecontainerchange', setup_form);
