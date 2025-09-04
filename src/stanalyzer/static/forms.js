/**
 * Given a string that looks like ``some[text][with][brackets]``, return the
 * part inside the nth set of square brackets `[]`
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
 * Select elements that include a given class
 *
 * @returns {jQuery} elements
 */
function jq_class(cls, rel="~=") {
    return $(sel_class(cls, rel));
}

/**
 * Return an class selector string
 *
 * @returns {String} CSS selector
 */
function sel_class(sel, rel='~=') {
    return sel_attr("class", sel, rel);
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
                case 'checkbox':
                    ref_value = ref_elem.prop('checked');
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
        if (['button', 'number', 'text'].includes(elem_type))
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
        if (toggle[0].checked) {
            target.show();
            target.filter('textarea').pushStack(
                target.find('textarea')).trigger('keyup');
        } else {
            target.hide();
        }
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

function result_span(message, result_type) {
    var error_template;
    var success_template;
    if (error_template === undefined)
        error_template = $('<span class="errmsg"></span>');
    if (success_template === undefined)
        success_template = $('<span class="success"></span>');

    let template = (result_type) ? success_template : error_template;
    const new_elem = template.clone();
    new_elem.text(message);
    return new_elem;
}

function error_span(message) {
    return result_span(message, false);
}

function success_span(message) {
    return result_span(message, true);
}

function clear_messages() {
    function rm_elem(idx, elem) {
        $(elem).remove();
    }
    jq_class('errmsg').each(rm_elem);
    jq_class('success').each(rm_elem);
    jq_class("sta-invalid").each((idx, elem) => $(elem).removeClass("sta-invalid"));
}

function delete_analysis(elem) {
    const form = $(`form#${active_page}_form`);
    const analysis_id = $(elem).data('id');
    form.attr('action', `/analysis/${analysis_id}`);
    form.data('id', analysis_id);
    submit_async(elem);
}

function rm_result(n_removed) {
    const form = $(`form#${active_page}_form`);
    const analysis_id = form.data('id');
    const to_remove = jq_id(`analysis_${analysis_id}`);
    n_removed = Number.parseInt(n_removed);

    // invalid int --> NaN, which fails all comparisons
    // thus, if n_removed > 0 is true, n_removed is a valid integer
    if (n_removed > 0) {
        to_remove.remove();
        update_global_message('<p>Job removed from database.</p><p style=\"font-size: .8em;\">You can still access your files on disk.</p>', true);
        return;
    }
    update_global_message("Failed to remove. An unexpected error has occurred.", false);
}

function button_outer(elem) {
    // returns a button's outermost element
    const parents = $(elem).parents(".ui-btn");
    if (!parents.length)
        return null;  // elem probably isn't a button

    return parents.last();
}

function popup_placeholder(elem) {
    const popup = $(elem).closest("[data-role=popup]");
    if (!popup.length)
        return null;  // elem probably isn't a popup

    const id = popup.attr('id');
    const placeholder = $(`#${id}-placeholder`);
    if (!placeholder.length)
        return null;

    return placeholder.prevAll("[data-rel=popup]");
}

function insert_before_original(elem, to_insert) {
    // Inserts something before whatever button started the interaction
    const insert_loc = popup_placeholder(elem) || button_outer(elem);
    if (insert_loc == null) {
        console.log("Can't find insert location");
        console.log(elem);
        return false;
    }
    return insert_loc.before(to_insert);
}

function submit_async(elem=null) {
    const form = $(`form#${active_page}_form`);

    // clear previous invalid results
    // jq_class("sta-invalid").each(elem => elem.removeClass("sta-invalid"));
    jq_class("sta-invalid").each((idx, elem) => $(elem).removeClass("sta-invalid"));

    clear_messages();
    const settings = {
        url: form.attr('action'),
        method: form.attr('method'),
        data: JSON.stringify(form.serializeJSON()),
        contentType: 'application/json',
        success: function(response) {  // default action; overridden via eval_settings
            update_global_message('Success', true);
        },
        error: undefined,
        statusCode: {
            422: function(response) {
                const data = response.responseJSON;
                const detail = data['detail'];

                if (detail === undefined)
                    return;

                var messages = [];
                $.each(detail, function(idx, item) {
                    // last field is the one where the error occurred
                    // other fields are
                    const name = item['loc'].slice(-1);
                    if (!name)
                        return;

                    const label = jq_label(name);
                    if (label !== undefined)
                        label.addClass("sta-invalid");

                    label.append(error_span(item['msg']));
                });
                update_global_message('Request could not be completed due to above errors.', false);
            },
            500: function(response) {
                update_global_message('An internal server error has occurred.', false);
            }
        }
    }

    const eval_settings = ["data", "success", "error"];

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

function analysis_failed(response_obj) {
    const response_json = response_obj.responseJSON;
    const is_success = false;

    if (response_json.length == 0) {
        update_global_message("No analysis selected", is_success);
        return;
    }

    const successes = [];
    const fails = [];
    for (const job of response_json) {
        if (job.pid === -1)
            fails.push(job.args);
        else
            successes.push(job.args);
    }

    if (fails.length == 0) {
        update_global_message("An unknown error has occurred", is_success);
        return;
    }

    const lines = [];
    if (successes.length) {
        lines.push("<div class=\"success\">");
        lines.push("<p>Analysis started with command(s):</p><pre>");
        lines.push(...successes);
        lines.push("</pre></div>");
    }

    lines.push("<p>Failed to start:</p><pre>");
    lines.push(...fails);
    lines.push("</pre>");

    update_global_message(lines.join("\n"), is_success);
}

function show_cmds(response_json) {
    // response_json should be an Array
    const msg = "<p>Analysis started with command(s):</p><pre>" +
        response_json.map(x => `${x['args']}`).join('\n') +
        "</pre><p style=\"font-size: .8em;\">(You'll still have to check your job's <code>.out/.err</code> files to see if the job succeeded)</p>";
    update_global_message(msg, true);
}

function update_projects_menu(response_json) {
    clear_messages();

    const request_type = response_json.request_type;
    const new_menu = response_json.data;
    const settings_elem = jq_id('json');
    const menu_elem = jq_name('project');

    // update in-form JSON representation
    settings_elem.val(JSON.stringify(new_menu));

    // remove old menu
    menu_elem.children(':not(option[value=add_new])').remove();

    // build new menu
    $.each(new_menu, function(id, obj) {
        const text = obj.title;
        menu_elem.prepend(`<option value="${id}">${text}</option>`);
    });

    switch (request_type) {
        case 'post':
            const new_title = jq_id('title').val();
            const new_item = menu_elem.children().filter(
                (idx, child) => child.innerHTML == new_title);
            menu_elem.val(new_item.val());
        case 'delete':
            menu_elem.trigger("change");
            break;
    }

    update_global_message('Success', true);
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

// popups w/o an activation link need to be explicitly initialized
$(document).on("DOMContentLoaded", function() {
    const msg = $('#global_message');
    const msg_content = $('#global_message-content');
    const msg_control = $('#message_control');
    const msg_close = msg.find('button');

    msg.popup();

    // jqm's dynamic positioning (via top/left properties) interferes with
    // style.css's positioning (via bottom/right)
    msg.data('mobile-popup')._reposition = $.noop;

    // setup other popup behavior
    $('#global_message-screen').on('click', function() {
        msg.popup("close");
        msg_control.show();
    });
    msg_control.on('click', function() {
        msg_control.hide();
        msg.popup("open");
    });
    msg_close.on('click', function() {
        msg.popup("close");
    })
});

function update_global_message(text, is_success=true, open=true) {
    var msg_container, msg_elem, msg_control;
    if (msg_elem === undefined) {
        msg_container = $('#global_message');
        msg_elem = $('#global_message-content');
        msg_control = $('#message_control');
    }

    msg_elem.html(text);

    msg_elem.removeClass(['errmsg', 'success']);

    const cls = (is_success) ? 'success' : 'errmsg';
    msg_elem.addClass(cls);

    if (open) {
        msg_container.popup("open");
        msg_control.hide();
    }
}
